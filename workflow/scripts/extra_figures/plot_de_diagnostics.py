from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import configure_matplotlib

configure_matplotlib()


if "snakemake" in locals():
    degs_path = Path(snakemake.input["degs"])
    outputs = {name: Path(value) for name, value in snakemake.output.items()}
    label_top = int(snakemake.params["label_top_genes"])
else:
    degs_path = Path("results/dea/interaction.tsv")
    outputs = {
        "csv": Path("results/extra_figures/source_data/de_diagnostics.tsv"),
        "summary": Path("results/extra_figures/de_significance_summary.tsv"),
        "pvalue_png": Path("plots/extra_figures/de_pvalue_histogram.png"),
        "pvalue_pdf": Path("plots/extra_figures/de_pvalue_histogram.pdf"),
        "ma_png": Path("plots/extra_figures/de_ma.png"),
        "ma_pdf": Path("plots/extra_figures/de_ma.pdf"),
        "volcano_png": Path("plots/extra_figures/de_volcano.png"),
        "volcano_pdf": Path("plots/extra_figures/de_volcano.pdf"),
        "counts_png": Path("plots/extra_figures/de_significant_counts.png"),
        "counts_pdf": Path("plots/extra_figures/de_significant_counts.pdf"),
    }
    label_top = 20

for path in outputs.values():
    path.parent.mkdir(parents=True, exist_ok=True)

degs = pd.read_csv(degs_path, sep="\t")
degs["minus_log10_pvalue"] = -np.log10(degs["pvalue"].clip(lower=np.finfo(float).tiny))
degs["label"] = degs["genesymbol"].fillna(degs["gene_id"])
degs.to_csv(outputs["csv"], sep="\t", index=False)

summary = pd.DataFrame(
    {
        "threshold": ["raw_pvalue_0.05", "bh_padj_0.05"],
        "genes": [int(degs["pvalue"].lt(0.05).sum()), int(degs["padj"].lt(0.05).sum())],
    }
)
summary.to_csv(outputs["summary"], sep="\t", index=False)

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.hist(degs["pvalue"].dropna(), bins=np.linspace(0, 1, 21), color="#4c78a8", edgecolor="white")
ax.set(xlabel="Raw p-value", ylabel="Genes", title="Interaction p-value distribution")
fig.tight_layout()
fig.savefig(outputs["pvalue_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["pvalue_pdf"], bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 5))
sns.scatterplot(
    data=degs, x="baseMean", y="log2FoldChange", hue="significant_padj_0.05",
    palette={False: "#808080", True: "#d62728"}, s=12, alpha=0.55, linewidth=0, ax=ax,
)
ax.set_xscale("log")
ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
ax.set(xlabel="Mean normalized count", ylabel="Interaction log2 fold-change", title="MA plot")
fig.tight_layout()
fig.savefig(outputs["ma_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["ma_pdf"], bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(degs["log2FoldChange"], degs["minus_log10_pvalue"], s=10, alpha=0.45, color="#4c78a8")
labels = degs.nsmallest(label_top, "pvalue").dropna(subset=["log2FoldChange", "minus_log10_pvalue"])
for row in labels.itertuples():
    ax.annotate(row.label, (row.log2FoldChange, row.minus_log10_pvalue), fontsize=8)
ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
ax.set(xlabel="Interaction log2 fold-change", ylabel="-log10(raw p-value)", title="Interaction volcano")
fig.tight_layout()
fig.savefig(outputs["volcano_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["volcano_pdf"], bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 4))
sns.barplot(data=summary, x="threshold", y="genes", color="#4c78a8", ax=ax)
for patch, value in zip(ax.patches, summary["genes"]):
    ax.text(patch.get_x() + patch.get_width() / 2, value, str(value), ha="center", va="bottom")
ax.set(xlabel="", ylabel="Genes", title="Significant interaction genes")
ax.tick_params(axis="x", rotation=20)
fig.tight_layout()
fig.savefig(outputs["counts_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["counts_pdf"], bbox_inches="tight")
plt.close(fig)
