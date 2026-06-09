# %% [markdown]
# # <center> Differential expression analysis
# %% [markdown]
# ## <center> Description
# %% [markdown]
# - env: squidpy_v06
# %% [markdown]
# ## <center> Set-up
# %% [markdown]
# #### **libs**
# %%
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import sys
from pathlib import Path
from scipy import sparse
## for vst
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# When this notebook is exported and run as a script, add the project root to sys.path
# so local workflow.* imports behave like they do in the notebook kernel.
if "__file__" in globals():
    PROJECT_ROOT = Path(__file__).resolve().parents[3]
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(PROJECT_ROOT))

## helper functions
from workflow.scripts.funs.helper_funs import get_project_root

## plotting
# import marsilea as ma
# import marsilea.plotter as mp
# from matplotlib.colors import NoNorm # for marsilea dotplot
import plotnine as p9
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# %% [markdown]
# #### **paths**
# %%
from pathlib import Path

if "snakemake" in globals():
    # in_metadata = Path(snakemake.input["in_metadata"])
    # out_metadata = Path(snakemake.output["out_metadata"])
    print('Hello World')
else:
    ## get paths
    ROOT = get_project_root(globals())
    ## input
    in_adata = ROOT / "output/preprocessing/objects/QC/adata.h5ad"
    ## output
    ### objects
    out_deg = ROOT / "output/DEA/objects/pyDESeq2/deg.tsv"
    out_deg.parent.mkdir(parents=True, exist_ok=True) # make sure the output folder exists
    ### plots
    out_plots = ROOT / "output/DEA/plots/pyDESeq2"
    out_plots.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# #### **reproducibility**
# %%
import random
import numpy as np
import torch
import scanpy.external as sce

SEED = 0

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)

# %% [markdown]
# ## <center> Main
# %% [markdown]
# #### **Load data**
# %%
adata = sc.read_h5ad(in_adata)

adata

# %%
# remove group B

adata = adata[adata.obs['QC_group']=='A'].copy()

adata

# %%
# make sure counts are raw

adata.X = adata.layers['cts.raw'].copy()

# %%
# check dist of samples across groups

adata.obs[['disease_group', 'condition_group']].value_counts()

# %%
# make sure the factors and levels are in the right order (disease_group and condition_group)

# These are the only levels used in the interaction model below.
# The first category is the reference level for each factor.
condition_levels = ["IFNG_neg", "IFNG_pos"]
disease_levels = ["HR_control", "HR_case"]

# Keep only samples belonging to the requested 2 x 2 design.
# This protects the interaction model from additional levels that are not part of the question.
adata = adata[
    adata.obs["condition_group"].isin(condition_levels)
    & adata.obs["disease_group"].isin(disease_levels)
].copy()

# Store the variables as categorical with explicit ordering.
# The formula below also encodes the reference levels explicitly for pyDESeq2/formulaic.
adata.obs["condition_group"] = pd.Categorical(
    adata.obs["condition_group"],
    categories=condition_levels,
)
adata.obs["disease_group"] = pd.Categorical(
    adata.obs["disease_group"],
    categories=disease_levels,
)

# Combined label used only for checking sample balance and plotting.
adata.obs["dea_group"] = pd.Categorical(
    adata.obs["condition_group"].astype(str) + "_" + adata.obs["disease_group"].astype(str),
    categories=[
        "IFNG_neg_HR_control",
        "IFNG_pos_HR_control",
        "IFNG_neg_HR_case",
        "IFNG_pos_HR_case",
    ],
)

# Sanity checks before model fitting.
assert adata.n_obs > 0, "No samples left after subsetting to the requested DEA groups."
assert not adata.obs[["condition_group", "disease_group"]].isna().any().any(), "Unexpected factor levels found."

dea_sample_counts = (
    adata.obs[["condition_group", "disease_group", "dea_group"]]
    .value_counts(sort=False)
    .rename("n_samples")
    .reset_index()
)

# A full-rank interaction model needs samples in each of the four cells.
assert (dea_sample_counts["n_samples"] > 0).all(), "The 2 x 2 interaction design has an empty group."

dea_sample_counts

# %% [markdown]
# ## <center> DEA
# %%
# prepare data for pyDESeq2

# Use the raw count layer prepared in QC.ipynb.
# Genes were already filtered there, so no additional expression filtering is applied here.
adata_dds = adata.copy()
adata_dds.X = adata_dds.layers["cts.raw"].copy()

# pyDESeq2 expects an in-memory dense integer count matrix.
# The rounding is a safety step in case the count layer was stored as floating point.
if sparse.issparse(adata_dds.X):
    adata_dds.X = adata_dds.X.toarray()
adata_dds.X = np.rint(np.asarray(adata_dds.X)).astype(int)

# Keep all QC-filtered genes, but remove any accidental all-zero genes because DESeq2 cannot estimate them.
nonzero_genes = np.asarray(adata_dds.X.sum(axis=0)).ravel() > 0
print(f"Genes before all-zero safety check: {adata_dds.n_vars}")
print(f"Genes after all-zero safety check: {nonzero_genes.sum()}")
adata_dds = adata_dds[:, nonzero_genes].copy()

adata_dds

# %%
# fit the DESeq2 model with the requested interaction design

# The C(..., Treatment(...)) terms encode the desired reference levels explicitly:
# - condition_group reference: IFNG_neg
# - disease_group reference: HR_control
# The interaction term lets the IFNG effect differ between HR_case and HR_control.
design_formula = (
    "~ C(condition_group, Treatment('IFNG_neg'))"
    " * C(disease_group, Treatment('HR_control'))"
)

dds = DeseqDataSet(
    adata=adata_dds,
    design=design_formula,
    refit_cooks=True,
    n_cpus=1,  # kept at 1 because this environment has shown pyDESeq2 multiprocessing issues
)

# Run the standard DESeq2 workflow: size factors, dispersions, LFCs, and Cook's refitting if needed.
dds.deseq2()

# Check the design matrix columns so the contrast below is easy to audit.
dds.obsm["design_matrix"].head()

# %%
# extract the requested contrast

# This contrast compares the fitted mean for IFNG_pos / HR_case against
# the fitted mean for IFNG_neg / HR_control. Therefore, a positive log2FC means
# higher expression in IFNG_pos HR_case relative to IFNG_neg HR_control.
contrast_ifng_hr = (
    dds.cond(condition_group="IFNG_pos", disease_group="HR_case")
    - dds.cond(condition_group="IFNG_neg", disease_group="HR_control")
).to_numpy()

stat_res = DeseqStats(
    dds,
    contrast=contrast_ifng_hr,
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
    n_cpus=1,
)
stat_res.summary()

# Format results and add gene symbols when they are available in adata.var.
deg = stat_res.results_df.copy()
deg.index.name = "gene_id"

if "genesymbol" in adata_dds.var.columns:
    deg.insert(0, "genesymbol", adata_dds.var["genesymbol"].reindex(deg.index).astype(str))

deg["contrast"] = "IFNG_pos_HR_case_vs_IFNG_neg_HR_control"
deg["significant_pvalue_0.05"] = deg["pvalue"] < 0.05
deg["significant_padj_0.05"] = deg["padj"] < 0.05
deg["direction"] = np.select(
    [
        deg["log2FoldChange"] > 0,
        deg["log2FoldChange"] < 0,
    ],
    ["higher_IFNG_pos_HR_case", "higher_IFNG_neg_HR_control"],
    default="no_change",
)

# Keep the most significant genes at the top of the saved table.
deg = deg.sort_values(["padj", "pvalue"], na_position="last")

deg.head(20)

# %%
# QC plot: P value histogram

# A broadly well-calibrated test often has a roughly uniform null component plus
# enrichment near zero when there is signal. Strong deviations can indicate model or batch issues.
plot_tab = deg.dropna(subset=["pvalue"]).copy()

w, h = 5, 4
p = (
    p9.ggplot(plot_tab, p9.aes(x="pvalue"))
    + p9.geom_histogram(
        breaks=np.arange(0, 1.05, 0.05),
        fill="#4C78A8",
        color="white",
        boundary=0,
    )
    + p9.scale_x_continuous(breaks=np.arange(0, 1.05, 0.1), limits=(0, 1))
    + p9.theme_classic()
    + p9.labs(
        title="P value distribution",
        x="Raw P value",
        y="Number of genes",
    )
    + p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=13),
        axis_title=p9.element_text(size=12),
        axis_text=p9.element_text(size=10),
        figure_size=(w, h),
    )
)

p.save(out_plots / "pvalue_histogram.png", width=w, height=h, limitsize=False)
p.save(out_plots / "pvalue_histogram.pdf", width=w, height=h, limitsize=False)

plt.close()

# %%
# QC plot: volcano plot with decoupler

# decoupler labels the top significant genes using the dataframe index.
# The y-axis uses adjusted P values so the significance threshold is FDR < 0.05.
# Use gene symbols for labels when available; otherwise fall back to gene IDs.
volcano_tab = deg.dropna(subset=["log2FoldChange", "padj"]).copy()
if "genesymbol" in volcano_tab.columns:
    volcano_tab["plot_label"] = volcano_tab["genesymbol"].where(
        volcano_tab["genesymbol"].notna() & (volcano_tab["genesymbol"] != "nan"),
        volcano_tab.index.astype(str),
    )
else:
    volcano_tab["plot_label"] = volcano_tab.index.astype(str)
volcano_tab = volcano_tab.set_index("plot_label", drop=False)

fig = dc.pl.volcano(
    data=volcano_tab,
    x="log2FoldChange",
    y="padj",
    top=20,
    thr_stat=0.5,
    thr_sign=0.05,
    figsize=(7, 6),
    dpi=120,
    return_fig=True,
)
fig.axes[0].set_title("IFNG_pos HR_case vs IFNG_neg HR_control")
fig.savefig(out_plots / "volcano.png", bbox_inches="tight")
fig.savefig(out_plots / "volcano.pdf", bbox_inches="tight")

plt.close(fig)

# %%
# QC plot: number of genes passing raw and adjusted P value thresholds

# Count genes passing the two requested significance thresholds.
sig_counts = pd.DataFrame(
    {
        "threshold": ["Raw P < 0.05", "Padj < 0.05"],
        "n_genes": [
            int((deg["pvalue"] < 0.05).sum()),
            int((deg["padj"] < 0.05).sum()),
        ],
    }
)

w, h = 4.5, 4
p = (
    p9.ggplot(sig_counts, p9.aes(x="threshold", y="n_genes", fill="threshold"))
    + p9.geom_col(color="black", width=0.7, show_legend=False)
    + p9.geom_text(p9.aes(label="n_genes"), va="bottom", size=10, format_string="{:.0f}")
    + p9.scale_fill_manual(["#4C78A8", "#F58518"])
    + p9.theme_classic()
    + p9.labs(
        title="Significant genes",
        x="",
        y="Number of genes",
    )
    + p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=13),
        axis_title=p9.element_text(size=12),
        axis_text_x=p9.element_text(size=10, rotation=20, ha="right"),
        axis_text_y=p9.element_text(size=10),
        figure_size=(w, h),
    )
)

p.save(out_plots / "significant_gene_counts.png", width=w, height=h, limitsize=False)
p.save(out_plots / "significant_gene_counts.pdf", width=w, height=h, limitsize=False)

plt.close()

sig_counts

# %%
# QC plot: MA plot

# The MA plot shows the relationship between average normalized expression and log2FC.
# It is useful for spotting low-count instability or global shifts.
ma_tab = deg.dropna(subset=["baseMean", "log2FoldChange"]).copy()
ma_tab["is_padj_0.05"] = ma_tab["padj"] < 0.05

w, h = 5.5, 4.5
p = (
    p9.ggplot(ma_tab, p9.aes(x="baseMean", y="log2FoldChange"))
    + p9.geom_point(p9.aes(color="is_padj_0.05"), alpha=0.55, size=0.8)
    + p9.scale_x_log10()
    + p9.scale_color_manual({False: "#7F7F7F", True: "#D62728"})
    + p9.geom_hline(yintercept=0, linetype="dashed")
    + p9.theme_classic()
    + p9.labs(
        title="MA plot",
        x="Mean normalized counts",
        y="log2FC",
        color="Padj < 0.05",
    )
    + p9.theme(
        plot_title=p9.element_text(hjust=0.5, size=13),
        axis_title=p9.element_text(size=12),
        axis_text=p9.element_text(size=10),
        legend_title=p9.element_text(size=10),
        legend_text=p9.element_text(size=9),
        figure_size=(w, h),
    )
)

p.save(out_plots / "ma_plot.png", width=w, height=h, limitsize=False)
p.save(out_plots / "ma_plot.pdf", width=w, height=h, limitsize=False)

plt.close()

# %% [markdown]
# #### **Save results**
# %%
deg.to_csv(out_deg,sep='\t')
