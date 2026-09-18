#!/usr/bin/env python

from pathlib import Path
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_wrap,
    geom_col,
    geom_errorbar,
    geom_jitter,
    ggplot,
    guide_legend,
    guides,
    labs,
    scale_fill_manual,
    theme,
    theme_classic,
)
from scipy import sparse


ROOT = Path(".")
IN_ADATA = ROOT / "output/preprocessing/objects/QC/adata.h5ad"
OUT_DIR = ROOT / "output/DEA/plots/selected_genes_expression_4groups_barpoints"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SELECTED_GENES = [
    "APOL1", "APOL6", "IRF1", "CCL2", "CXCL9", "CXCL10",
    "ICAM1", "IDO1", "PDCD1LG2", "GBP1", "PSMB9", "CCDC68",
    "CTSS", "TMEM140", "TRIM21", "STAT1", "STAT2", "STAT3",
    "IRF2", "IRF8", "IRF9", "CXCL11",
]

GROUP_LABELS = {
    "HR_case__IFNG_neg": "HR case, IFNG-",
    "HR_case__IFNG_pos": "HR case, IFNG+",
    "HR_control__IFNG_neg": "HR control, IFNG-",
    "HR_control__IFNG_pos": "HR control, IFNG+",
}

GROUP_ORDER = [
    "HR case, IFNG-",
    "HR case, IFNG+",
    "HR control, IFNG-",
    "HR control, IFNG+",
]

GROUP_COLORS = {
    "HR case, IFNG-": "#fdb863",
    "HR case, IFNG+": "#e66101",
    "HR control, IFNG-": "#92c5de",
    "HR control, IFNG+": "#2166ac",
}


def make_unique(names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    unique = []
    for name in names:
        if name not in counts:
            counts[name] = 0
            unique.append(name)
        else:
            counts[name] += 1
            unique.append(f"{name}.{counts[name]}")
    return unique


def build_plot(plot_df: pd.DataFrame, summary_df: pd.DataFrame, error_col: str, title: str):
    plots_per_row = 8
    n_rows = int(np.ceil(len(summary_df["genesymbol"].cat.categories) / plots_per_row))
    width = 2.3 * plots_per_row
    height = 4 * n_rows

    summary_tmp = summary_df.copy()
    summary_tmp["ymin"] = summary_tmp["mean_expr"] - summary_tmp[error_col]
    summary_tmp["ymax"] = summary_tmp["mean_expr"] + summary_tmp[error_col]

    p = (
        ggplot(summary_tmp, aes(x="grp", y="mean_expr", fill="grp"))
        + geom_col(width=0.66, color="black", alpha=0.65, show_legend=False)
        + geom_errorbar(
            aes(ymin="ymin", ymax="ymax"),
            width=0.18,
            size=0.5,
            show_legend=False,
        )
        + geom_jitter(
            data=plot_df,
            mapping=aes(x="grp", y="expr", fill="grp"),
            color="black",
            size=2.3,
            width=0.10,
            alpha=0.95,
            show_legend=True,
        )
        + scale_fill_manual(values=GROUP_COLORS)
        + guides(fill=guide_legend(override_aes={"size": 5}))
        + theme_classic()
        + facet_wrap("~genesymbol", scales="free_y", ncol=plots_per_row)
        + labs(title=title, x=None, y="expr", fill="grp")
        + theme(
            plot_title=element_text(ha="center", size=14),
            strip_text=element_text(size=12),
            axis_title=element_text(size=14),
            axis_text_x=element_blank(),
            legend_text=element_text(size=11),
            legend_title=element_text(size=12),
            figure_size=(width, height),
        )
    )
    return p, width, height


adata = ad.read_h5ad(IN_ADATA)
adata = adata[adata.obs["QC_group"] == "A"].copy()

if "cts.deseq2_norm" not in adata.layers:
    raise ValueError("adata.h5ad does not contain the 'cts.deseq2_norm' layer.")

adata.obs = adata.obs.copy()
adata.var = adata.var.copy()

renamed = []
for ens, symbol in zip(adata.var_names, adata.var["genesymbol"]):
    if pd.isna(symbol) or symbol == "":
        renamed.append(str(ens))
    else:
        renamed.append(str(symbol))
adata.var_names = make_unique(renamed)

present_genes = [gene for gene in SELECTED_GENES if gene in adata.var_names]
missing_genes = [gene for gene in SELECTED_GENES if gene not in adata.var_names]
if missing_genes:
    warnings.warn(
        "Missing selected genes in normalized expression object: "
        + ", ".join(missing_genes)
    )

expr_mat = adata[:, present_genes].layers["cts.deseq2_norm"]
if sparse.issparse(expr_mat):
    expr_mat = expr_mat.toarray()
expr_df = pd.DataFrame(expr_mat, index=adata.obs_names, columns=present_genes)
expr_df.index.name = "sample_id"

plot_df = (
    expr_df.reset_index()
    .melt(id_vars="sample_id", var_name="genesymbol", value_name="expr")
    .merge(adata.obs[["grp"]], left_on="sample_id", right_index=True, how="left")
)
plot_df["grp"] = plot_df["grp"].map(GROUP_LABELS)
plot_df = plot_df.dropna(subset=["grp"]).copy()
plot_df["grp"] = pd.Categorical(plot_df["grp"], categories=GROUP_ORDER, ordered=True)
plot_df["genesymbol"] = pd.Categorical(plot_df["genesymbol"], categories=present_genes, ordered=True)

summary_df = (
    plot_df.groupby(["genesymbol", "grp"], observed=True)["expr"]
    .agg(["mean", "std", "count"])
    .reset_index()
    .rename(columns={"mean": "mean_expr", "std": "sd_expr", "count": "n"})
)
summary_df["sem_expr"] = summary_df["sd_expr"] / np.sqrt(summary_df["n"])
summary_df["genesymbol"] = pd.Categorical(summary_df["genesymbol"], categories=present_genes, ordered=True)
summary_df["grp"] = pd.Categorical(summary_df["grp"], categories=GROUP_ORDER, ordered=True)

p_sd, width, height = build_plot(
    plot_df,
    summary_df,
    error_col="sd_expr",
    title="Expression of selected genes (mean ± SD with all samples)",
)
p_sd.save(
    OUT_DIR / "selected_genes_expression_4groups_mean_sd.png",
    width=width,
    height=height,
    limitsize=False,
    dpi=300,
)
p_sd.save(
    OUT_DIR / "selected_genes_expression_4groups_mean_sd.pdf",
    width=width,
    height=height,
    limitsize=False,
)

p_sem, _, _ = build_plot(
    plot_df,
    summary_df,
    error_col="sem_expr",
    title="Expression of selected genes (mean ± SEM with all samples)",
)
p_sem.save(
    OUT_DIR / "selected_genes_expression_4groups_mean_sem.png",
    width=width,
    height=height,
    limitsize=False,
    dpi=300,
)
p_sem.save(
    OUT_DIR / "selected_genes_expression_4groups_mean_sem.pdf",
    width=width,
    height=height,
    limitsize=False,
)
