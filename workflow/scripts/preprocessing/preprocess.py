import random
from pathlib import Path

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from pydeseq2.dds import DeseqDataSet


if "snakemake" in locals():
    input_path = Path(snakemake.input["adata"])
    output_path = Path(snakemake.output["adata"])
    counts_path = Path(snakemake.output["counts"])
    summary_path = Path(snakemake.output["summary"])
    seed = int(snakemake.params["seed"])
    min_cells = int(snakemake.params["min_cells"])
    min_total_count = int(snakemake.params["min_total_count"])
    min_count = int(snakemake.params["min_count"])
    min_prop = float(snakemake.params["min_prop"])
    highly_variable_genes = int(snakemake.params["highly_variable_genes"])
    pca_components = int(snakemake.params["pca_components"])
    neighbors = int(snakemake.params["neighbors"])
    leiden_resolution = float(snakemake.params["leiden_resolution"])
else:
    input_path = Path("results/preprocessing/adata_unfiltered.h5ad")
    output_path = Path("results/preprocessing/adata.h5ad")
    counts_path = Path("results/preprocessing/counts.tsv")
    summary_path = Path("results/preprocessing/qc_summary.tsv")
    seed = 0
    min_cells, min_total_count, min_count, min_prop = 5, 10, 10, 0.25
    highly_variable_genes, pca_components, neighbors, leiden_resolution = 2000, 20, 15, 0.15

for path in [output_path, counts_path, summary_path]:
    path.parent.mkdir(parents=True, exist_ok=True)
random.seed(seed)
np.random.seed(seed)

adata = sc.read_h5ad(input_path)
adata.obs["grp"] = (
    adata.obs["disease_group"].astype(str) + "__" + adata.obs["condition_group"].astype(str)
)
genes_initial = adata.n_vars
sc.pp.filter_genes(adata, min_cells=min_cells)
genes_after_min_cells = adata.n_vars
dc.pp.filter_by_expr(
    adata=adata,
    min_total_count=min_total_count,
    group="grp",
    min_prop=min_prop,
    min_count=min_count,
    large_n=0,
    inplace=True,
)
genes_after_expression = adata.n_vars

adata.X.data = np.rint(adata.X.data)
adata.layers["cts.raw"] = adata.X.copy()
adata_dds = adata.copy()
adata_dds.X = adata_dds.X.toarray()
dds = DeseqDataSet(adata=adata_dds, design="~ 1", refit_cooks=True, n_cpus=1)
dds.vst(use_design=False)
adata.layers["cts.vst"] = np.asarray(dds.layers["vst_counts"])
adata.layers["cts.deseq2_norm"] = np.asarray(dds.layers["normed_counts"])
adata.X = adata.layers["cts.deseq2_norm"].copy()
sc.pp.log1p(adata)
adata.layers["cts.deseq2_norm_log"] = np.asarray(adata.X).copy()

sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=highly_variable_genes)
sc.tl.pca(
    adata,
    layer="cts.vst",
    n_comps=pca_components,
    use_highly_variable=True,
    svd_solver="arpack",
    random_state=seed,
)
sc.pp.neighbors(adata, n_neighbors=neighbors, n_pcs=pca_components, random_state=seed)
sc.tl.leiden(
    adata,
    key_added=f"leiden_{leiden_resolution}",
    resolution=leiden_resolution,
    random_state=seed,
    flavor="igraph",
    n_iterations=2,
)
sc.tl.umap(adata, random_state=seed)

if adata.n_obs != 66:
    raise ValueError(f"Expected 66 total samples, found {adata.n_obs}")
qc_a = adata.obs["QC_group"].astype(str).eq("A")
if int(qc_a.sum()) != 50:
    raise ValueError(f"Expected 50 QC-group A samples, found {int(qc_a.sum())}")
complete_pairs = (
    adata.obs.loc[qc_a]
    .groupby("individualID", observed=True)["condition_group"]
    .nunique()
    .eq(2)
)
if int(complete_pairs.sum()) != 25:
    raise ValueError(f"Expected 25 complete QC-group A donor pairs, found {int(complete_pairs.sum())}")

counts = pd.DataFrame(
    adata.layers["cts.raw"].toarray(), index=adata.obs_names, columns=adata.var_names
).T
counts.to_csv(counts_path, sep="\t", index_label="gene_id")
pd.DataFrame(
    {
        "metric": [
            "samples_total",
            "samples_qc_group_a",
            "complete_qc_group_a_pairs",
            "genes_initial",
            "genes_after_min_cells",
            "genes_after_expression_filter",
        ],
        "value": [
            adata.n_obs,
            int(qc_a.sum()),
            int(complete_pairs.sum()),
            genes_initial,
            genes_after_min_cells,
            genes_after_expression,
        ],
    }
).to_csv(summary_path, sep="\t", index=False)
adata.write_h5ad(output_path)
