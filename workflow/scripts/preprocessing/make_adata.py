import sys
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import read_checked_tsv


if "snakemake" in locals():
    metadata_path = Path(snakemake.input["metadata"])
    counts_path = Path(snakemake.input["counts"])
    annotation_path = Path(snakemake.input["network"])
    output_path = Path(snakemake.output["adata"])
    expected_sha256 = snakemake.params["sha256"]
else:
    metadata_path = Path("results/preprocessing/metadata.tsv")
    counts_path = Path("data/gene_counts.txt")
    annotation_path = Path("data/prior_knowledge/ensembl_gene_symbols_20260608.tsv")
    output_path = Path("results/preprocessing/adata_unfiltered.h5ad")
    expected_sha256 = "07278e787896b42a0487e835a1a88817b9405f234bd45ed425724e7df5a55d22"

output_path.parent.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(metadata_path, sep="\t").set_index("sample_id")
counts = pd.read_csv(counts_path, sep="\t").set_index("Gene").T
annotation = read_checked_tsv(
    annotation_path,
    expected_sha256,
    ["gene_id", "gene_symbol", "gene_biotype", "chromosome"],
).drop_duplicates("gene_id").set_index("gene_id")

missing_samples = metadata.index.difference(counts.index)
if len(missing_samples):
    raise ValueError(f"Count matrix is missing samples: {missing_samples.tolist()}")
counts = counts.reindex(metadata.index)

adata = ad.AnnData(
    X=sparse.csr_matrix(counts.astype("float32")),
    obs=metadata.copy(),
    var=pd.DataFrame(index=counts.columns),
)
adata.obs_names.name = "sample_id"
adata.var_names.name = "gene_id"
adata.var = adata.var.join(annotation, how="left")
adata.obs["disease_group"] = pd.Categorical(
    adata.obs["disease_group"], categories=["HR_control", "HR_case"], ordered=True
)
adata.obs["QC_group"] = pd.Categorical(
    adata.obs["QC_group"], categories=["A", "B"], ordered=True
)
adata.obs["condition_group"] = pd.Categorical(
    adata.obs["condition_group"], categories=["IFNG_neg", "IFNG_pos"], ordered=True
)
adata.obs["individualID"] = adata.obs["individualID"].astype("category")

if adata.n_obs != 66:
    raise ValueError(f"Expected 66 samples, found {adata.n_obs}")
adata.write_h5ad(output_path)
