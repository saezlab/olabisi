from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


if "snakemake" in locals():
    enrichment_path = Path(snakemake.input["csv"])
    collectri_csv = Path(snakemake.output["collectri_csv"])
    collectri_png = Path(snakemake.output["collectri_png"])
    collectri_pdf = Path(snakemake.output["collectri_pdf"])
    progeny_csv = Path(snakemake.output["progeny_csv"])
    progeny_png = Path(snakemake.output["progeny_png"])
    progeny_pdf = Path(snakemake.output["progeny_pdf"])
    top_n = int(snakemake.params["top_n"])
else:
    enrichment_path = Path("results/enrichment/ulm.tsv")
    collectri_csv = Path("results/paper/source_data/Fig5E_collectri.tsv")
    collectri_png = Path("plots/paper/Fig5E_collectri.png")
    collectri_pdf = Path("plots/paper/Fig5E_collectri.pdf")
    progeny_csv = Path("results/paper/source_data/Fig5F_progeny.tsv")
    progeny_png = Path("plots/paper/Fig5F_progeny.png")
    progeny_pdf = Path("plots/paper/Fig5F_progeny.pdf")
    top_n = 10

for path in [collectri_csv, collectri_png, collectri_pdf, progeny_csv, progeny_png, progeny_pdf]:
    path.parent.mkdir(parents=True, exist_ok=True)

enrichment = pd.read_csv(enrichment_path, sep="\t")


def plot_overview(gene_set, title, csv_path, png_path, pdf_path):
    table = (
        enrichment[enrichment["gene_set"].eq(gene_set)]
        .sort_values(["padj", "source"], kind="stable")
        .head(top_n)
        .copy()
    )
    table["neglog10_padj_x_direction"] = (
        -np.log10(table["padj"].clip(lower=np.finfo(float).tiny))
        * np.sign(table["activity"])
    )
    table = table.sort_values("activity", kind="stable")
    table.to_csv(csv_path, sep="\t", index=False)

    norm = mpl.colors.TwoSlopeNorm(
        vmin=min(float(table["activity"].min()), -1e-9),
        vcenter=0,
        vmax=max(float(table["activity"].max()), 1e-9),
    )
    cmap = mpl.colormaps["coolwarm"]
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.barh(
        table["source"],
        table["neglog10_padj_x_direction"],
        color=cmap(norm(table["activity"])),
        edgecolor="black",
        linewidth=0.7,
        alpha=0.8,
    )
    threshold = -np.log10(0.05)
    ax.axvline(0, color="black", linewidth=0.7)
    ax.axvline(threshold, color="black", linestyle="--", linewidth=0.7)
    ax.axvline(-threshold, color="black", linestyle="--", linewidth=0.7)
    ax.set_title(title)
    ax.set_xlabel("signed -log10(BH-adjusted p-value)")
    ax.set_ylabel("")
    scalar = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(scalar, ax=ax, label="ULM activity")
    fig.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


plot_overview("collectri", "Top CollecTRI TFs (ULM)", collectri_csv, collectri_png, collectri_pdf)
plot_overview("progeny", "Top PROGENy pathways (ULM)", progeny_csv, progeny_png, progeny_pdf)
