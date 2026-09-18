import random
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse


if "snakemake" in locals():
    adata_path = Path(snakemake.input["adata"])
    output_path = Path(snakemake.output["degs"])
    seed = int(snakemake.params["seed"])
else:
    adata_path = Path("results/preprocessing/adata.h5ad")
    output_path = Path("results/dea/interaction.tsv")
    seed = 0

output_path.parent.mkdir(parents=True, exist_ok=True)
random.seed(seed)
np.random.seed(seed)

adata = sc.read_h5ad(adata_path)
adata = adata[adata.obs["QC_group"].astype(str).eq("A")].copy()
adata.X = adata.layers["cts.raw"].copy()
adata.obs["ifng_pos"] = adata.obs["condition_group"].astype(str).eq("IFNG_pos").astype(int)
adata.obs["case_ifng_pos"] = (
    adata.obs["condition_group"].astype(str).eq("IFNG_pos")
    & adata.obs["disease_group"].astype(str).eq("HR_case")
).astype(int)

if adata.n_obs != 50:
    raise ValueError(f"Expected 50 QC-group A samples, found {adata.n_obs}")
pair_counts = adata.obs.groupby("individualID", observed=True)["condition_group"].nunique()
if len(pair_counts) != 25 or not pair_counts.eq(2).all():
    raise ValueError("Paired interaction analysis requires 25 donors with both IFNG conditions")

adata_dds = adata.copy()
adata_dds.X = adata_dds.layers["cts.raw"].copy()
if sparse.issparse(adata_dds.X):
    adata_dds.X = adata_dds.X.toarray()
adata_dds.X = np.rint(np.asarray(adata_dds.X)).astype(int)
nonzero = np.asarray(adata_dds.X.sum(axis=0)).ravel() > 0
adata_dds = adata_dds[:, nonzero].copy()
adata_dds.obs["individualID"] = adata_dds.obs["individualID"].astype("category")

design_formula = "~ individualID + ifng_pos + case_ifng_pos"
dds = DeseqDataSet(adata=adata_dds, design=design_formula, refit_cooks=True, n_cpus=1)
dds.deseq2()
design_columns = dds.obsm["design_matrix"].columns.tolist()
if "case_ifng_pos" not in design_columns:
    raise ValueError(f"Interaction coefficient missing from design matrix: {design_columns}")
contrast = np.zeros(len(design_columns), dtype=float)
contrast[design_columns.index("case_ifng_pos")] = 1.0
stats = DeseqStats(
    dds,
    contrast=contrast,
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True,
    n_cpus=1,
)
stats.summary()

degs = stats.results_df.copy()
degs.index.name = "gene_id"
degs = degs.reset_index()
degs.insert(1, "genesymbol", adata_dds.var["gene_symbol"].reindex(degs["gene_id"]).to_numpy())
degs["contrast"] = "interaction_(HR_case_IFNGpos_vs_IFNGneg)_vs_(HR_control_IFNGpos_vs_IFNGneg)"
degs["design_formula"] = design_formula
degs["coefficient_tested"] = "case_ifng_pos"
degs["significant_pvalue_0.05"] = degs["pvalue"].lt(0.05)
degs["significant_padj_0.05"] = degs["padj"].lt(0.05)
degs["direction"] = np.select(
    [degs["log2FoldChange"].gt(0), degs["log2FoldChange"].lt(0)],
    ["stronger_IFNG_response_in_HR_case", "stronger_IFNG_response_in_HR_control"],
    default="no_change",
)
degs = degs.sort_values(["padj", "pvalue", "gene_id"], na_position="last", kind="stable")
degs.to_csv(output_path, sep="\t", index=False)
