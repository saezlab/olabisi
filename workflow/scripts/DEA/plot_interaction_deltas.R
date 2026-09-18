suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

in_paired <- "output/DEA/objects/pyDESeq2_interaction_diagnostics/top20_gene_paired_deltas.tsv"
in_deg <- "output/DEA/objects/pyDESeq2_interaction/deg.tsv"
out_dir <- "output/DEA/plots/pyDESeq2_interaction_diagnostics"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paired <- read_tsv(in_paired, show_col_types = FALSE)
deg <- read_tsv(in_deg, show_col_types = FALSE)

top_genes <- deg %>%
  slice_min(order_by = pvalue, n = 12, with_ties = FALSE) %>%
  select(gene_id, genesymbol, pvalue)

plot_df <- paired %>%
  inner_join(top_genes, by = "gene_id") %>%
  mutate(
    gene_label = if_else(!is.na(genesymbol) & genesymbol != "nan", genesymbol, gene_id),
    gene_label = reorder(gene_label, pvalue),
    disease_group = factor(disease_group, levels = c("HR_control", "HR_case"))
  )

p <- ggplot(plot_df, aes(x = gene_label, y = delta_ifng, color = disease_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 2) +
  stat_summary(
    aes(group = disease_group),
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.5),
    size = 3,
    shape = 18
  ) +
  facet_wrap(~ disease_group, ncol = 1) +
  scale_color_manual(values = c("HR_control" = "#1f78b4", "HR_case" = "#d95f02")) +
  labs(
    title = "Donor-level IFNG response deltas for top interaction genes",
    subtitle = "Each point is IFNG_pos - IFNG_neg within one donor; diamonds mark group medians",
    x = NULL,
    y = "Delta raw counts",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "donor_ifng_deltas_top12.png"), p, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "donor_ifng_deltas_top12.pdf"), p, width = 10, height = 7)
