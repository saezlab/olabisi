suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

case_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_case__IFNGpos_vs_IFNGneg/deg.tsv"
control_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_control__IFNGpos_vs_IFNGneg/deg.tsv"
interaction_path <- "output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"
out_dir <- "output/DEA/plots/log2fc_scatter_ifng_contrasts"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    genesymbol = coalesce(genesymbol_control, genesymbol_case, genesymbol_interaction)
  )

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
  scale_color_gradient2(
    low = "#1f78b4",
    mid = "grey92",
    high = "#d95f02",
    midpoint = 0,
    limits = c(-color_limit, color_limit),
    name = "Interaction\nlog2FC"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits) +
  labs(
    title = "log2FC comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- log2FC",
    y = "HR case IFNG+ vs IFNG- log2FC",
    caption = paste(
      "Point color shows interaction-term log2FC; zero genes showed significantly different IFNG response",
      "in case vs controls (disease-treatment interaction, padj < 0.05).",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  file.path(out_dir, "ifng_contrast_log2fc_scatter_by_interaction.png"),
  p,
  width = 8,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_log2fc_scatter_by_interaction.pdf"),
  p,
  width = 8,
  height = 8
)
