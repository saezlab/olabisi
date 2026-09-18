#!/usr/bin/env python

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.image import imread


ROOT = Path(".")
out_pdf = ROOT / "output/DEA/reports/interaction_analysis_report.pdf"
out_pdf.parent.mkdir(parents=True, exist_ok=True)

deg = pd.read_csv(ROOT / "output/DEA/objects/pyDESeq2_interaction/deg.tsv", sep="\t")
cal_summary = pd.read_csv(
    ROOT / "output/DEA/objects/calibrate_pvalues_interaction/calibration_summary.tsv",
    sep="\t",
)

top_hits = deg.loc[:, ["gene_id", "genesymbol", "log2FoldChange", "pvalue", "padj", "direction"]].head(10).copy()
top_hits["log2FoldChange"] = top_hits["log2FoldChange"].map(lambda x: f"{x:.3f}")
top_hits["pvalue"] = top_hits["pvalue"].map(lambda x: f"{x:.3g}")
top_hits["padj"] = top_hits["padj"].map(lambda x: f"{x:.3g}")

summary_map = dict(zip(cal_summary["metric"], cal_summary["value"]))


def add_text_page(pdf, title, lines, fontsize=11, title_size=18):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.patch.set_facecolor("white")
    plt.axis("off")
    left = 0.11
    fig.text(left, 0.94, title, fontsize=title_size, fontweight="bold", va="top")
    y = 0.885
    for block in lines:
        wrapped = textwrap.fill(block, width=95) if block else ""
        fig.text(left, y, wrapped, fontsize=fontsize, va="top", family="monospace" if block.startswith("~ ") else None)
        y -= 0.045 + 0.013 * wrapped.count("\n")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def add_image_page(pdf, title, image_paths, captions):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.patch.set_facecolor("white")
    left = 0.11
    fig.text(left, 0.94, title, fontsize=18, fontweight="bold", va="top")
    positions = [(0.11, 0.55, 0.78, 0.24), (0.11, 0.16, 0.78, 0.24)]
    for (img_path, caption, pos) in zip(image_paths, captions, positions):
        ax = fig.add_axes(pos)
        ax.imshow(imread(img_path))
        ax.axis("off")
        fig.text(pos[0], pos[1] + pos[3] + 0.01, caption, fontsize=11, va="bottom")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


with PdfPages(out_pdf) as pdf:
    add_text_page(
        pdf,
        "Interaction Analysis Report",
        [
            "Objective: test whether the IFNG response differs between HR_case and HR_control in the paired PBMC dataset.",
            "Input object: output/preprocessing/objects/QC/adata.h5ad",
            "Analysis outputs: output/DEA/objects/pyDESeq2_interaction and output/DEA/objects/calibrate_pvalues_interaction",
            "Sample structure after QC group A restriction: 25 paired donors, 50 samples total.",
            "Group balance: HR_control IFNG_neg = 13, HR_control IFNG_pos = 13, HR_case IFNG_neg = 12, HR_case IFNG_pos = 12.",
        ],
    )

    add_text_page(
        pdf,
        "Model Specification",
        [
            "Three formulations were relevant in this analysis.",
            "~ DEA_group",
            "Original one-factor model in the existing PyDESeq2 script. This compares combined 4-level groups directly, for example HR_case__IFNG_pos vs HR_control__IFNG_pos. It does not directly test a contrast-of-contrasts.",
            "~ individualID + condition_group + disease_group:condition_group",
            "Initial paired interaction model that directly encodes the scientific question. This is the standard contrast-of-contrasts specification, but it produced a singular design in this PyDESeq2 workflow on the current dataset.",
            "Briefly: disease_group is fixed within donor, so the donor indicators already encode case/control status. The categorical interaction expansion added one redundant case-by-condition column, making the matrix rank-deficient.",
            "~ individualID + ifng_pos + case_ifng_pos",
            "Final fitted parameterization. Here ifng_pos is the within-donor IFNG effect in the reference setting and case_ifng_pos is the additional IFNG effect in cases. The coefficient case_ifng_pos is the tested interaction term.",
        ],
    )

    add_text_page(
        pdf,
        "Interaction Results Summary",
        [
            f"Design fitted: {deg['design_formula'].iloc[0]}",
            f"Coefficient tested: {deg['coefficient_tested'].iloc[0]}",
            f"Genes tested: {len(deg)}",
            f"Raw p < 0.05: {int((deg['pvalue'] < 0.05).sum())}",
            f"Raw padj < 0.05: {int((deg['padj'] < 0.05).sum())}",
            "Top interaction hits by raw p-value:",
            top_hits.to_string(index=False),
        ],
        fontsize=10,
    )

    add_image_page(
        pdf,
        "Primary Diagnostic Plots",
        [
            ROOT / "output/DEA/plots/pyDESeq2_interaction/volcano.png",
            ROOT / "output/DEA/plots/pyDESeq2_interaction/ma_plot.png",
        ],
        [
            "Volcano plot for the fitted interaction coefficient.",
            "MA plot for the fitted interaction coefficient.",
        ],
    )

    add_image_page(
        pdf,
        "P-value Distributions",
        [
            ROOT / "output/DEA/plots/pyDESeq2_interaction/pvalue_histogram.png",
            ROOT / "output/DEA/plots/calibrate_pvalues_interaction/P_hist_calibrated.png",
        ],
        [
            "Raw interaction-model p-value histogram.",
            "Calibrated p-value histogram using the repo's fdrtool-based calibration workflow.",
        ],
    )

    add_text_page(
        pdf,
        "P-value Calibration",
        [
            "Calibration followed the same approach already used elsewhere in the repository: fdrtool applied to the Wald statistics assuming a normal null.",
            f"Estimated null SD from fdrtool: {summary_map['null_sd']:.4f}",
            f"Raw p < 0.05: {int(summary_map['raw_p_lt_0.05'])}",
            f"Raw padj < 0.05: {int(summary_map['raw_padj_lt_0.05'])}",
            f"Calibrated p < 0.05: {int(summary_map['cal_p_lt_0.05'])}",
            f"Calibrated padj < 0.05: {int(summary_map['cal_padj_lt_0.05'])}",
            "Interpretation: the null SD below 1 indicates conservative statistics. Calibration increases the number of nominally significant genes, but no genes pass FDR < 0.05 after calibration.",
        ],
    )
