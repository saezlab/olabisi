from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import configure_matplotlib

configure_matplotlib()


if "snakemake" in locals():
    adata_path = Path(snakemake.input["adata"])
    heatmap_csv = Path(snakemake.output["heatmap_csv"])
    heatmap_png = Path(snakemake.output["heatmap_png"])
    heatmap_pdf = Path(snakemake.output["heatmap_pdf"])
    umap_csv = Path(snakemake.output["umap_csv"])
    groups_png = Path(snakemake.output["groups_png"])
    groups_pdf = Path(snakemake.output["groups_pdf"])
    markers_png = Path(snakemake.output["markers_png"])
    markers_pdf = Path(snakemake.output["markers_pdf"])
    podocyte_markers = list(snakemake.params["podocyte_markers"])
    umap_markers = list(snakemake.params["umap_markers"])
else:
    adata_path = Path("results/preprocessing/adata.h5ad")
    heatmap_csv = Path("results/paper/source_data/Fig5B_heatmap.tsv")
    heatmap_png = Path("plots/paper/Fig5B_heatmap.png")
    heatmap_pdf = Path("plots/paper/Fig5B_heatmap.pdf")
    umap_csv = Path("results/paper/source_data/Fig5C_umap.tsv")
    groups_png = Path("plots/paper/Fig5C_umap_groups.png")
    groups_pdf = Path("plots/paper/Fig5C_umap_groups.pdf")
    markers_png = Path("plots/paper/Fig5C_umap_markers.png")
    markers_pdf = Path("plots/paper/Fig5C_umap_markers.pdf")
    podocyte_markers = [
        "APOL1", "ACTN4", "LAMB2", "ITGB1", "SYNPO", "VEGFA", "KIRREL1",
        "COL4A5", "ITGA3", "LAMA5", "SULF1", "NPNT", "CD2AP", "TCF21",
        "FAT1", "PODXL", "PLA2R1", "MAFB",
    ]
    umap_markers = ["PLA2R1", "ITGA3", "PODXL", "OLFM3"]

for path in [heatmap_csv, heatmap_png, heatmap_pdf, umap_csv, groups_png, groups_pdf, markers_png, markers_pdf]:
    path.parent.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(adata_path)
symbols = adata.var["gene_symbol"].astype("string")
symbol_to_index = {}
for index, symbol in zip(adata.var_names, symbols):
    if pd.notna(symbol) and symbol not in symbol_to_index:
        symbol_to_index[str(symbol)] = index

missing_heatmap = sorted(set(podocyte_markers).difference(symbol_to_index))
missing_umap = sorted(set(umap_markers).difference(symbol_to_index))
if missing_heatmap or missing_umap:
    raise ValueError(
        f"Missing configured markers; heatmap={missing_heatmap}, umap={missing_umap}"
    )

log_values = np.asarray(adata.layers["cts.deseq2_norm_log"])
gene_position = {gene_id: position for position, gene_id in enumerate(adata.var_names)}
heatmap_matrix = np.column_stack(
    [log_values[:, gene_position[symbol_to_index[symbol]]] for symbol in podocyte_markers]
)
heatmap_source = pd.DataFrame(heatmap_matrix, index=adata.obs_names, columns=podocyte_markers)
heatmap_source.index.name = "sample_id"
heatmap_long = heatmap_source.reset_index().melt(
    id_vars="sample_id", var_name="genesymbol", value_name="log1p_gene_expression"
)
heatmap_long = heatmap_long.merge(
    adata.obs[["condition_group", "disease_group", "QC_group"]].reset_index(),
    on="sample_id",
    how="left",
    validate="many_to_one",
)
heatmap_long.to_csv(heatmap_csv, sep="\t", index=False)

fig, ax = plt.subplots(figsize=(10, 4.8))
sns.heatmap(
    heatmap_source.T,
    cmap="magma",
    linewidths=0.1,
    linecolor="black",
    cbar_kws={"label": "gene expression"},
    ax=ax,
)
ax.set_title("Expression of selected podocyte markers across samples")
ax.set_xlabel("Sample")
ax.set_ylabel("")
ax.tick_params(axis="x", labelsize=8)
ax.tick_params(axis="y", labelsize=8)
fig.tight_layout()
fig.savefig(heatmap_png, dpi=300, bbox_inches="tight")
fig.savefig(heatmap_pdf, bbox_inches="tight")
plt.close(fig)

umap_source = pd.DataFrame(
    adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]
).join(adata.obs[["individualID", "condition_group", "disease_group", "QC_group"]])
for symbol in umap_markers:
    umap_source[symbol] = log_values[:, gene_position[symbol_to_index[symbol]]]
umap_source.index.name = "sample_id"
umap_source.reset_index().to_csv(umap_csv, sep="\t", index=False)

condition_colors = {"IFNG_neg": "#1f77b4", "IFNG_pos": "#ff7f0e"}
disease_markers = {"HR_control": "o", "HR_case": "*"}
fig, ax = plt.subplots(figsize=(6, 5))
for disease, marker in disease_markers.items():
    for condition, color in condition_colors.items():
        mask = (
            umap_source["disease_group"].astype(str).eq(disease)
            & umap_source["condition_group"].astype(str).eq(condition)
        )
        ax.scatter(
            umap_source.loc[mask, "UMAP1"],
            umap_source.loc[mask, "UMAP2"],
            s=70,
            marker=marker,
            facecolor=color,
            edgecolor="black",
            linewidth=0.4,
        )
legend = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="black", label=condition)
    for condition, color in condition_colors.items()
] + [
    Line2D([0], [0], marker=marker, color="black", linestyle="none", label=disease)
    for disease, marker in disease_markers.items()
]
ax.legend(handles=legend, frameon=False, fontsize=8)
ax.set_title("Unsupervised clustering of samples")
ax.set_xlabel("UMAP1")
ax.set_ylabel("UMAP2")
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(groups_png, dpi=300, bbox_inches="tight")
fig.savefig(groups_pdf, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
for ax, symbol in zip(axes.flat, umap_markers):
    points = ax.scatter(
        umap_source["UMAP1"],
        umap_source["UMAP2"],
        c=umap_source[symbol],
        cmap="viridis",
        s=45,
        edgecolor="black",
        linewidth=0.25,
    )
    ax.set_title(symbol)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    fig.colorbar(points, ax=ax, label="log1p normalized expression")
fig.savefig(markers_png, dpi=300, bbox_inches="tight")
fig.savefig(markers_pdf, bbox_inches="tight")
plt.close(fig)
