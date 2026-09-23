from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import configure_matplotlib

configure_matplotlib()


if "snakemake" in locals():
    adata_path = Path(snakemake.input["adata"])
    outputs = {name: Path(value) for name, value in snakemake.output.items()}
    n_components = int(snakemake.params["pca_components"])
else:
    adata_path = Path("results/preprocessing/adata.h5ad")
    outputs = {
        "pca_csv": Path("results/extra_figures/source_data/pca.tsv"),
        "anova": Path("results/extra_figures/source_data/pca_anova.tsv"),
        "umap_csv": Path("results/extra_figures/source_data/umap_qc.tsv"),
        "pca_png": Path("plots/extra_figures/pca_qc.png"),
        "pca_pdf": Path("plots/extra_figures/pca_qc.pdf"),
        "umap_png": Path("plots/extra_figures/umap_qc.png"),
        "umap_pdf": Path("plots/extra_figures/umap_qc.pdf"),
    }
    n_components = 6

for path in outputs.values():
    path.parent.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(adata_path)
metadata_columns = ["individualID", "disease_group", "condition_group", "QC_group", "leiden_0.15"]
metadata = adata.obs[metadata_columns].astype(str).reset_index(names="sample_id")

pca = pd.DataFrame(
    np.asarray(adata.obsm["X_pca"])[:, :n_components],
    columns=[f"PC{i + 1}" for i in range(n_components)],
)
pca = pd.concat([metadata.reset_index(drop=True), pca], axis=1)
pca.to_csv(outputs["pca_csv"], sep="\t", index=False)

anova_rows = []
for component in [f"PC{i + 1}" for i in range(n_components)]:
    for factor in ["disease_group", "condition_group", "QC_group", "individualID"]:
        groups = [
            values[component].to_numpy()
            for _, values in pca.groupby(factor, observed=True)
            if len(values) > 1
        ]
        if len(groups) < 2:
            continue
        statistic, pvalue = stats.f_oneway(*groups)
        grand_mean = pca[component].mean()
        between = sum(len(group) * (group.mean() - grand_mean) ** 2 for group in groups)
        total = ((pca[component] - grand_mean) ** 2).sum()
        anova_rows.append(
            {
                "component": component,
                "factor": factor,
                "groups": len(groups),
                "f_statistic": statistic,
                "pvalue": pvalue,
                "eta_squared": between / total if total else np.nan,
            }
        )
pd.DataFrame(anova_rows).sort_values(["component", "pvalue", "factor"]).to_csv(
    outputs["anova"], sep="\t", index=False
)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, factor in zip(axes, ["disease_group", "condition_group"]):
    sns.scatterplot(data=pca, x="PC1", y="PC2", hue=factor, style="QC_group", s=65, ax=ax)
    ax.set_title(f"PCA by {factor.replace('_', ' ')}")
    ax.axhline(0, color="0.85", linewidth=0.7)
    ax.axvline(0, color="0.85", linewidth=0.7)
fig.tight_layout()
fig.savefig(outputs["pca_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["pca_pdf"], bbox_inches="tight")
plt.close(fig)

umap = pd.DataFrame(np.asarray(adata.obsm["X_umap"]), columns=["UMAP1", "UMAP2"])
umap = pd.concat([metadata.reset_index(drop=True), umap], axis=1)
umap.to_csv(outputs["umap_csv"], sep="\t", index=False)

fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, factor in zip(
    axes.flat, ["disease_group", "condition_group", "QC_group", "leiden_0.15"]
):
    sns.scatterplot(data=umap, x="UMAP1", y="UMAP2", hue=factor, s=60, ax=ax)
    ax.set_title(factor.replace("_", " "))
    ax.set_xticks([])
    ax.set_yticks([])
fig.tight_layout()
fig.savefig(outputs["umap_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["umap_pdf"], bbox_inches="tight")
plt.close(fig)
