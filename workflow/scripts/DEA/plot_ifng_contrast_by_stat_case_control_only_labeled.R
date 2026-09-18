suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

case_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_case__IFNGpos_vs_IFNGneg/deg.tsv"
control_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_control__IFNGpos_vs_IFNGneg/deg.tsv"
interaction_path <- "output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"
out_dir <- "output/DEA/plots/log2fc_scatter_ifng_contrasts"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jak_stat_genes <- c(
  "APOL1", "APOL6", "IRF1", "CCL2", "CXCL9", "CXCL10",
  "ICAM1", "IDO1", "PDCD1LG2", "GBP1", "PSMB9", "CCDC68",
  "CTSS", "TMEM140", "TRIM21", "STAT1", "STAT2", "STAT3",
  "IRF2", "IRF8", "IRF9", "CXCL11"
)

read_deg <- function(path, prefix) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(gene_id, genesymbol, stat, log2FoldChange, padj) %>%
    rename(
      !!paste0("genesymbol_", prefix) := genesymbol,
      !!paste0("stat_", prefix) := stat,
      !!paste0("log2FoldChange_", prefix) := log2FoldChange,
      !!paste0("padj_", prefix) := padj
    )
}

case_deg <- read_deg(case_path, "case")
control_deg <- read_deg(control_path, "control")
interaction_deg <- read_tsv(interaction_path, show_col_types = FALSE) %>%
  select(gene_id, genesymbol, stat, log2FoldChange) %>%
  rename(
    genesymbol_interaction = genesymbol,
    stat_interaction = stat,
    log2FoldChange_interaction = log2FoldChange
  )

plot_df <- control_deg %>%
  inner_join(case_deg, by = "gene_id") %>%
  inner_join(interaction_deg, by = "gene_id") %>%
  mutate(
    genesymbol = coalesce(genesymbol_control, genesymbol_case, genesymbol_interaction),
    label_gene = !is.na(genesymbol) & genesymbol %in% jak_stat_genes,
    case_hit = !is.na(padj_case) & padj_case < 0.05 & abs(log2FoldChange_case) > 2,
    control_hit = !is.na(padj_control) & padj_control < 0.05 & abs(log2FoldChange_control) > 2,
    case_only = case_hit & !control_hit,
    control_only = control_hit & !case_hit
  )

label_df <- plot_df %>%
  filter(label_gene)

missing_labels <- setdiff(jak_stat_genes, label_df$genesymbol)
if (length(missing_labels) > 0) {
  warning("Selected genes not found in joined data: ", paste(missing_labels, collapse = ", "))
}

axis_limits <- range(c(plot_df$stat_control, plot_df$stat_case), na.rm = TRUE)

make_plot <- function(highlight_col, highlight_color, highlight_label, with_labels = TRUE) {
  highlight_flag <- plot_df[[highlight_col]]
  background_df <- plot_df[!highlight_flag, , drop = FALSE]
  highlight_df <- plot_df[highlight_flag, , drop = FALSE]

  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_point(
      data = background_df,
      aes(x = stat_control, y = stat_case),
      color = "grey75",
      alpha = 0.8,
      size = 1.6
    ) +
    geom_point(
      data = highlight_df,
      aes(x = stat_control, y = stat_case),
      color = highlight_color,
      alpha = 0.9,
      size = 1.8
    ) +
    coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
    labs(
      title = "Stat comparison for IFNG+ vs IFNG- by disease group",
      subtitle = paste0(highlight_label, " genes: padj < 0.05 and |log2FC| > 2 in ", tolower(highlight_label), " but not the other contrast"),
      x = "HR control IFNG+ vs IFNG- stat",
      y = "HR case IFNG+ vs IFNG- stat",
      caption = if (with_labels) "Labels mark JAK/STAT-responsive genes." else NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0),
      plot.margin = margin(10, 35, 10, 10)
    )

  if (with_labels) {
    p <- p +
      geom_text_repel(
        data = label_df,
        aes(x = stat_control, y = stat_case, label = genesymbol),
        size = 3.2,
        color = "black",
        min.segment.length = 0,
        box.padding = 0.3,
        point.padding = 0.15,
        seed = 42,
        max.overlaps = Inf
      )
  }

  p
}

make_case_only_interaction_plot <- function(value_col, legend_title, subtitle_text) {
  highlight_df <- plot_df[plot_df$case_only, , drop = FALSE]
  color_limit <- max(abs(plot_df[[value_col]]), na.rm = TRUE)

  ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_point(
      data = highlight_df,
      aes(x = stat_control, y = stat_case, color = .data[[value_col]]),
      alpha = 0.9,
      size = 1.8
    ) +
    scale_color_gradient2(
      low = "#1f78b4",
      mid = "grey92",
      high = "#d95f02",
      midpoint = 0,
      limits = c(-color_limit, color_limit),
      name = legend_title
    ) +
    coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
    labs(
      title = "Stat comparison for IFNG+ vs IFNG- by disease group",
      subtitle = subtitle_text,
      x = "HR control IFNG+ vs IFNG- stat",
      y = "HR case IFNG+ vs IFNG- stat"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.margin = margin(10, 35, 10, 10)
    )
}

case_plot <- make_plot("case_only", "#d95f02", "Case only")
control_plot <- make_plot("control_only", "#1f78b4", "Control only")
case_plot_unlabeled <- make_plot("case_only", "#d95f02", "Case only", with_labels = FALSE)
control_plot_unlabeled <- make_plot("control_only", "#1f78b4", "Control only", with_labels = FALSE)
case_plot_interaction_stat <- make_case_only_interaction_plot(
  "stat_interaction",
  "Interaction\nstat",
  "Case-only genes only: color shows interaction stat"
)
case_plot_interaction_log2fc <- make_case_only_interaction_plot(
  "log2FoldChange_interaction",
  "Interaction\nlog2FC",
  "Case-only genes: color shows interaction log2FC; grey points are all other genes"
)

ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_labeled.png"),
  case_plot,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_labeled.pdf"),
  case_plot,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_control_only_hits_labeled.png"),
  control_plot,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_control_only_hits_labeled.pdf"),
  control_plot,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits.png"),
  case_plot_unlabeled,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits.pdf"),
  case_plot_unlabeled,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_control_only_hits.png"),
  control_plot_unlabeled,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_control_only_hits.pdf"),
  control_plot_unlabeled,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_by_interaction_stat.png"),
  case_plot_interaction_stat,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_by_interaction_stat.pdf"),
  case_plot_interaction_stat,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_by_interaction_log2fc.png"),
  case_plot_interaction_log2fc,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_case_only_hits_by_interaction_log2fc.pdf"),
  case_plot_interaction_log2fc,
  width = 9,
  height = 8
)
