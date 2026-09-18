# %% [markdown]
# # <center> Differential expression analysis: interaction model
# %% [markdown]
# Tests whether the IFNG response differs between HR_case and HR_control.
# With paired samples, this is implemented as:
#     ~ individualID + condition_group + disease_group:condition_group
# The interaction coefficient is the contrast-of-contrasts:
#     (HR_case IFNG_pos - HR_case IFNG_neg) - (HR_control IFNG_pos - HR_control IFNG_neg)

# %%
import random
import sys
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import scanpy as sc
from scipy import sparse

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

if "__file__" in globals():
    PROJECT_ROOT = Path(__file__).resolve().parents[3]
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(PROJECT_ROOT))

from workflow.scripts.funs.helper_funs import get_project_root


if "snakemake" in globals():
    in_adata = Path(snakemake.input["in_adata"])
    out_deg_interaction = Path(snakemake.output["out_deg_interaction"])
    out_deg_ifng_pos = Path(snakemake.output["out_deg_ifng_pos"])
    for out_deg in [out_deg_interaction, out_deg_ifng_pos]:
        out_deg.parent.mkdir(parents=True, exist_ok=True)
    out_plots_interaction = Path(snakemake.params["out_plots_interaction"])
    out_plots_ifng_pos = Path(snakemake.params["out_plots_ifng_pos"])
    out_plots_interaction.mkdir(parents=True, exist_ok=True)
    out_plots_ifng_pos.mkdir(parents=True, exist_ok=True)
else:
    ROOT = get_project_root(globals())
    in_adata = ROOT / "output/preprocessing/objects/QC/adata.h5ad"
    out_deg_interaction = ROOT / "output/DEA/objects/pyDESeq2_interaction/deg.tsv"
    out_deg_ifng_pos = ROOT / "output/DEA/objects/pyDESeq2_ifng_pos/deg.tsv"
    for out_deg in [out_deg_interaction, out_deg_ifng_pos]:
        out_deg.parent.mkdir(parents=True, exist_ok=True)
    out_plots_interaction = ROOT / "output/DEA/plots/pyDESeq2_interaction"
    out_plots_ifng_pos = ROOT / "output/DEA/plots/pyDESeq2_ifng_pos"
    out_plots_interaction.mkdir(parents=True, exist_ok=True)
    out_plots_ifng_pos.mkdir(parents=True, exist_ok=True)


SEED = 0
random.seed(SEED)
np.random.seed(SEED)

INDIVIDUAL_ID_MAP = {
    "SD1": "individual_01", "SD2": "individual_01", "SD3": "individual_02", "SD4": "individual_02",
    "SD5": "individual_03", "SD6": "individual_03", "SD7": "individual_04", "SD8": "individual_04",
    "SD9": "individual_05", "SD10": "individual_05", "SD11": "individual_06", "SD12": "individual_06",
    "SD13": "individual_07", "SD14": "individual_07", "SD15": "individual_08", "SD16": "individual_08",
    "SD17": "individual_09", "SD18": "individual_09", "SD19": "individual_10", "SD20": "individual_10",
    "SD21": "individual_11", "SD22": "individual_11", "SD23": "individual_12", "SD24": "individual_12",
    "SD25": "individual_13", "SD26": "individual_13", "SD27": "individual_14", "SD28": "individual_14",
    "SD29": "individual_15", "SD30": "individual_15", "SD31": "individual_16", "SD32": "individual_16",
    "SD33": "individual_17", "SD34": "individual_17", "SD35": "individual_18", "SD36": "individual_18",
    "SD37": "individual_19", "SD38": "individual_19", "SD39": "individual_20", "SD40": "individual_20",
    "SD41": "individual_21", "SD42": "individual_21", "SD43": "individual_22", "SD44": "individual_22",
    "SD45": "individual_23", "SD46": "individual_23", "SD47": "individual_24", "SD48": "individual_24",
    "SD49": "individual_25", "SD50": "individual_25", "SD51": "individual_26", "SD52": "individual_26",
    "SD53": "individual_27", "SD54": "individual_27", "SD55": "individual_28", "SD56": "individual_28",
    "SD57": "individual_29", "SD58": "individual_29", "SD59": "individual_30", "SD60": "individual_30",
    "SD61": "individual_31", "SD62": "individual_31", "SD63": "individual_32", "SD64": "individual_32",
    "SD65": "individual_33", "SD66": "individual_33",
}


adata = sc.read_h5ad(in_adata)

# Match the existing DEA preprocessing.
adata = adata[adata.obs["QC_group"] == "A"].copy()
adata.X = adata.layers["cts.raw"].copy()
adata.obs["individualID"] = adata.obs_names.map(INDIVIDUAL_ID_MAP)

condition_levels = ["IFNG_neg", "IFNG_pos"]
disease_levels = ["HR_control", "HR_case"]

adata = adata[
    adata.obs["condition_group"].isin(condition_levels)
    & adata.obs["disease_group"].isin(disease_levels)
].copy()

adata.obs["condition_group"] = pd.Categorical(
    adata.obs["condition_group"],
    categories=condition_levels,
)
adata.obs["disease_group"] = pd.Categorical(
    adata.obs["disease_group"],
    categories=disease_levels,
)

assert adata.n_obs > 0, "No samples left after subsetting to the requested 2 x 2 design."
assert not adata.obs[["condition_group", "disease_group"]].isna().any().any(), "Unexpected factor levels found."
assert adata.obs["individualID"].notna().all(), "Some samples could not be mapped to individualID for paired analysis."

# Encode the interaction explicitly to avoid redundant Patsy columns in the paired design.
adata.obs["ifng_pos"] = (adata.obs["condition_group"] == "IFNG_pos").astype(int)
adata.obs["case_ifng_pos"] = (
    (adata.obs["condition_group"] == "IFNG_pos")
    & (adata.obs["disease_group"] == "HR_case")
).astype(int)

paired_counts = (
    adata.obs[["individualID", "condition_group"]]
    .value_counts(sort=False)
    .rename("n_samples")
    .reset_index()
)
complete_pairs = paired_counts.groupby("individualID")["condition_group"].nunique().eq(len(condition_levels))
assert complete_pairs.all(), "Every individualID must have both IFNG_neg and IFNG_pos samples for the paired interaction model."

adata_dds = adata.copy()
adata_dds.X = adata_dds.layers["cts.raw"].copy()

if sparse.issparse(adata_dds.X):
    adata_dds.X = adata_dds.X.toarray()
adata_dds.X = np.rint(np.asarray(adata_dds.X)).astype(int)

nonzero_genes = np.asarray(adata_dds.X.sum(axis=0)).ravel() > 0
print(f"Genes before all-zero safety check: {adata_dds.n_vars}")
print(f"Genes after all-zero safety check: {nonzero_genes.sum()}")
adata_dds = adata_dds[:, nonzero_genes].copy()

adata_dds.obs["individualID"] = adata_dds.obs["individualID"].astype("category")
design_formula = "~ individualID + ifng_pos + case_ifng_pos"
coefficient_configs = {
    "case_ifng_pos": {
        "out_deg": out_deg_interaction,
        "out_plots": out_plots_interaction,
        "contrast_label": "interaction_(HR_case_IFNGpos_vs_IFNGneg)_vs_(HR_control_IFNGpos_vs_IFNGneg)",
        "pos_direction": "stronger_IFNG_response_in_HR_case",
        "neg_direction": "stronger_IFNG_response_in_HR_control",
        "volcano_title": "IFNG interaction: HR_case vs HR_control",
        "ma_ylabel": "Interaction log2FC",
    },
    "ifng_pos": {
        "out_deg": out_deg_ifng_pos,
        "out_plots": out_plots_ifng_pos,
        "contrast_label": "main_effect_IFNGpos_vs_IFNGneg_in_HR_control_reference",
        "pos_direction": "higher_in_IFNG_pos",
        "neg_direction": "higher_in_IFNG_neg",
        "volcano_title": "Main IFNG effect: IFNG_pos vs IFNG_neg",
        "ma_ylabel": "IFNG log2FC",
    },
}

print(f"Using design: {design_formula}")

dds = DeseqDataSet(
    adata=adata_dds,
    design=design_formula,
    refit_cooks=True,
    n_cpus=1,
)
dds.deseq2()
design_columns = dds.obsm["design_matrix"].columns.tolist()
print(design_columns)

def extract_and_plot(coefficient_name, config):
    contrast_vector = np.zeros(len(design_columns), dtype=float)
    contrast_vector[design_columns.index(coefficient_name)] = 1.0

    stat_res = DeseqStats(
        dds,
        contrast=contrast_vector,
        alpha=0.05,
        cooks_filter=True,
        independent_filter=True,
        n_cpus=1,
    )
    stat_res.summary()

    deg = stat_res.results_df.copy()
    deg.index.name = "gene_id"

    if "genesymbol" in adata_dds.var.columns:
        deg.insert(0, "genesymbol", adata_dds.var["genesymbol"].reindex(deg.index).astype(str))

    deg["contrast"] = config["contrast_label"]
    deg["design_formula"] = design_formula
    deg["coefficient_tested"] = coefficient_name
    deg["significant_pvalue_0.05"] = deg["pvalue"] < 0.05
    deg["significant_padj_0.05"] = deg["padj"] < 0.05
    deg["direction"] = np.select(
        [
            deg["log2FoldChange"] > 0,
            deg["log2FoldChange"] < 0,
        ],
        [config["pos_direction"], config["neg_direction"]],
        default="no_change",
    )

    deg = deg.sort_values(["padj", "pvalue"], na_position="last")

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
        + p9.labs(title="P value distribution", x="Raw P value", y="Number of genes")
        + p9.theme(
            plot_title=p9.element_text(hjust=0.5, size=13),
            axis_title=p9.element_text(size=12),
            axis_text=p9.element_text(size=10),
            figure_size=(w, h),
        )
    )
    p.save(config["out_plots"] / "pvalue_histogram.png", width=w, height=h, limitsize=False)
    p.save(config["out_plots"] / "pvalue_histogram.pdf", width=w, height=h, limitsize=False)
    plt.close()

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
    fig.axes[0].set_title(config["volcano_title"])
    fig.savefig(config["out_plots"] / "volcano.png", bbox_inches="tight")
    fig.savefig(config["out_plots"] / "volcano.pdf", bbox_inches="tight")
    plt.close(fig)

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
        + p9.labs(title="Significant genes", x="", y="Number of genes")
        + p9.theme(
            plot_title=p9.element_text(hjust=0.5, size=13),
            axis_title=p9.element_text(size=12),
            axis_text_x=p9.element_text(size=10, rotation=20, ha="right"),
            axis_text_y=p9.element_text(size=10),
            figure_size=(w, h),
        )
    )
    p.save(config["out_plots"] / "significant_gene_counts.png", width=w, height=h, limitsize=False)
    p.save(config["out_plots"] / "significant_gene_counts.pdf", width=w, height=h, limitsize=False)
    plt.close()

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
            y=config["ma_ylabel"],
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
    p.save(config["out_plots"] / "ma_plot.png", width=w, height=h, limitsize=False)
    p.save(config["out_plots"] / "ma_plot.pdf", width=w, height=h, limitsize=False)
    plt.close()

    deg.to_csv(config["out_deg"], sep="\t")


for coefficient_name, config in coefficient_configs.items():
    extract_and_plot(coefficient_name, config)
