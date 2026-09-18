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
    select(gene_id, genesymbol, log2FoldChange) %>%
    rename(
      !!paste0("genesymbol_", prefix) := genesymbol,
      !!paste0("log2FoldChange_", prefix) := log2FoldChange
    )
}

case_deg <- read_deg(case_path, "case")
control_deg <- read_deg(control_path, "control")
interaction_deg <- read_tsv(interaction_path, show_col_types = FALSE) %>%
  select(gene_id, genesymbol, log2FoldChange) %>%
  rename(
    genesymbol_interaction = genesymbol,
    log2FoldChange_interaction = log2FoldChange
  )

plot_df <- control_deg %>%
  inner_join(case_deg, by = "gene_id") %>%
  inner_join(interaction_deg, by = "gene_id") %>%
  mutate(
    genesymbol = coalesce(genesymbol_control, genesymbol_case, genesymbol_interaction),
    label_gene = !is.na(genesymbol) & genesymbol %in% jak_stat_genes
  )

label_df <- plot_df %>%
  filter(label_gene)

missing_labels <- setdiff(jak_stat_genes, label_df$genesymbol)
if (length(missing_labels) > 0) {
  warning("Selected genes not found in joined data: ", paste(missing_labels, collapse = ", "))
}

axis_limits <- range(
  c(plot_df$log2FoldChange_control, plot_df$log2FoldChange_case),
  na.rm = TRUE
)

color_limit <- max(abs(plot_df$log2FoldChange_interaction), na.rm = TRUE)

p <- ggplot(
  plot_df,
  aes(
    x = log2FoldChange_control,
    y = log2FoldChange_case,
    color = log2FoldChange_interaction
  )
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_point(alpha = 0.85, size = 1.6) +
  geom_text_repel(
    data = label_df,
    aes(label = genesymbol),
    size = 3.2,
    color = "black",
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.15,
    seed = 42,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_gradient2(
    low = "#1f78b4",
    mid = "grey92",
    high = "#d95f02",
    midpoint = 0,
    limits = c(-color_limit, color_limit),
    name = "Interaction\nlog2FC"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "log2FC comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- log2FC",
    y = "HR case IFNG+ vs IFNG- log2FC",
    caption = paste(
      "Point color shows interaction-term log2FC; zero genes showed significantly different IFNG response",
      "in case vs controls (disease-treatment interaction, padj < 0.05). Labels mark JAK/STAT-responsive genes.",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(10, 35, 10, 10)
  )

ggsave(
  file.path(out_dir, "ifng_contrast_log2fc_scatter_by_interaction_labeled.png"),
  p,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_log2fc_scatter_by_interaction_labeled.pdf"),
  p,
  width = 9,
  height = 8
)
