suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
})

case_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_case__IFNGpos_vs_IFNGneg/deg.tsv"
control_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_control__IFNGpos_vs_IFNGneg/deg.tsv"
interaction_path <- "output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"
ifng_pos_path <- "output/DEA/objects/calibrate_pvalues_ifng_pos/deg_calibrated.tsv"
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
    select(gene_id, genesymbol, stat) %>%
    rename(
      !!paste0("genesymbol_", prefix) := genesymbol,
      !!paste0("stat_", prefix) := stat
    )
}

case_deg <- read_deg(case_path, "case")
control_deg <- read_deg(control_path, "control")
interaction_deg <- read_tsv(interaction_path, show_col_types = FALSE) %>%
  select(gene_id, genesymbol, stat) %>%
  rename(
    genesymbol_interaction = genesymbol,
    stat_interaction = stat
  )
ifng_pos_deg <- read_tsv(ifng_pos_path, show_col_types = FALSE) %>%
  select(gene_id, genesymbol, stat) %>%
  rename(
    genesymbol_ifng_pos = genesymbol,
    stat_ifng_pos = stat
  )

plot_df <- control_deg %>%
  inner_join(case_deg, by = "gene_id") %>%
  inner_join(interaction_deg, by = "gene_id") %>%
  inner_join(ifng_pos_deg, by = "gene_id") %>%
  mutate(
    genesymbol = coalesce(
      genesymbol_control,
      genesymbol_case,
      genesymbol_interaction,
      genesymbol_ifng_pos
    ),
    label_gene = !is.na(genesymbol) & genesymbol %in% jak_stat_genes
  )

label_df <- plot_df %>%
  filter(label_gene)

missing_labels <- setdiff(jak_stat_genes, label_df$genesymbol)
if (length(missing_labels) > 0) {
  warning("Selected genes not found in joined data: ", paste(missing_labels, collapse = ", "))
}

axis_limits <- range(
  c(plot_df$stat_control, plot_df$stat_case),
  na.rm = TRUE
)

color_limit <- max(abs(plot_df$stat_interaction), na.rm = TRUE)
color_limit_ifng_pos <- max(abs(plot_df$stat_ifng_pos), na.rm = TRUE)

p <- ggplot(
  plot_df,
  aes(
    x = stat_control,
    y = stat_case,
    color = stat_interaction
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
    name = "Interaction\nstat"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "Stat comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- stat",
    y = "HR case IFNG+ vs IFNG- stat",
    caption = paste(
      "Point color shows interaction-term stat; zero genes showed significantly different IFNG response",
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

p_unlabeled <- ggplot(
  plot_df,
  aes(
    x = stat_control,
    y = stat_case,
    color = stat_interaction
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
    name = "Interaction\nstat"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "Stat comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- stat",
    y = "HR case IFNG+ vs IFNG- stat",
    caption = paste(
      "Point color shows interaction-term stat; zero genes showed significantly different IFNG response",
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

p_ifng_pos_unlabeled <- ggplot(
  plot_df,
  aes(
    x = stat_control,
    y = stat_case,
    color = stat_ifng_pos
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
    limits = c(-color_limit_ifng_pos, color_limit_ifng_pos),
    name = "ifng_pos\nstat"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "Stat comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- stat",
    y = "HR case IFNG+ vs IFNG- stat",
    caption = "Point color shows the ifng_pos-term stat from the same interaction model."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(10, 35, 10, 10)
  )

p_unlabeled_shared_ifng_pos_scale <- ggplot(
  plot_df,
  aes(
    x = stat_control,
    y = stat_case,
    color = stat_interaction
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
    limits = c(-color_limit_ifng_pos, color_limit_ifng_pos),
    name = "Interaction\nstat"
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits, clip = "off") +
  labs(
    title = "Stat comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- stat",
    y = "HR case IFNG+ vs IFNG- stat",
    caption = "Point color shows interaction-term stat; color scale matched to the ifng_pos-term stat plot for direct visual comparison."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(hjust = 0),
    plot.margin = margin(10, 35, 10, 10)
  )

ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_labeled.png"),
  p,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_labeled.pdf"),
  p,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat.png"),
  p_unlabeled,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat.pdf"),
  p_unlabeled,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_by_ifng_pos_stat.png"),
  p_ifng_pos_unlabeled,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_by_ifng_pos_stat.pdf"),
  p_ifng_pos_unlabeled,
  width = 9,
  height = 8
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_shared_ifng_pos_scale.png"),
  p_unlabeled_shared_ifng_pos_scale,
  width = 9,
  height = 8,
  dpi = 300
)
ggsave(
  file.path(out_dir, "ifng_contrast_by_stat_shared_ifng_pos_scale.pdf"),
  p_unlabeled_shared_ifng_pos_scale,
  width = 9,
  height = 8
)
