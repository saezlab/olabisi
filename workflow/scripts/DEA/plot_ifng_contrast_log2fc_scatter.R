suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

case_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_case__IFNGpos_vs_IFNGneg/deg.tsv"
control_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_control__IFNGpos_vs_IFNGneg/deg.tsv"
out_dir <- "output/DEA/plots/log2fc_scatter_ifng_contrasts"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_deg <- function(path, prefix) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(gene_id, genesymbol, log2FoldChange, padj) %>%
    rename(
      !!paste0("genesymbol_", prefix) := genesymbol,
      !!paste0("log2FoldChange_", prefix) := log2FoldChange,
      !!paste0("padj_", prefix) := padj
    )
}

case_deg <- read_deg(case_path, "case")
control_deg <- read_deg(control_path, "control")

plot_df <- control_deg %>%
  inner_join(case_deg, by = "gene_id") %>%
  mutate(
    genesymbol = coalesce(genesymbol_control, genesymbol_case),
    sig_control = !is.na(padj_control) & padj_control < 0.05 & abs(log2FoldChange_control) > 2,
    sig_case = !is.na(padj_case) & padj_case < 0.05 & abs(log2FoldChange_case) > 2,
    sig_group = case_when(
      sig_control & sig_case ~ "Both",
      sig_control ~ "Control only",
      sig_case ~ "Case only",
      TRUE ~ "Neither"
    ),
    sig_group = factor(sig_group, levels = c("Neither", "Control only", "Case only", "Both"))
  )

axis_limits <- range(
  c(plot_df$log2FoldChange_control, plot_df$log2FoldChange_case),
  na.rm = TRUE
)

p_plain <- ggplot(
  plot_df,
  aes(x = log2FoldChange_control, y = log2FoldChange_case)
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_point(alpha = 0.8, size = 1.6, color = "grey55") +
  coord_equal(xlim = axis_limits, ylim = axis_limits) +
  labs(
    title = "log2FC comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- log2FC",
    y = "HR case IFNG+ vs IFNG- log2FC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank()
  )

both_n <- sum(plot_df$sig_group == "Both")
case_only_n <- sum(plot_df$sig_group == "Case only")
control_only_n <- sum(plot_df$sig_group == "Control only")

p_colored <- ggplot(
  plot_df,
  aes(x = log2FoldChange_control, y = log2FoldChange_case, color = sig_group)
) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_point(alpha = 0.8, size = 1.6) +
  scale_color_manual(
    values = c(
      "Neither" = "grey75",
      "Control only" = "#1f78b4",
      "Case only" = "#d95f02",
      "Both" = "#6a3d9a"
    ),
    drop = FALSE
  ) +
  coord_equal(xlim = axis_limits, ylim = axis_limits) +
  labs(
    title = "log2FC comparison for IFNG+ vs IFNG- by disease group",
    x = "HR control IFNG+ vs IFNG- log2FC",
    y = "HR case IFNG+ vs IFNG- log2FC",
    color = "padj < 0.05\nand |log2FC| > 2",
    caption = paste0(
      "Both: ", both_n,
      "; Case only: ", case_only_n,
      "; Control only: ", control_only_n, "."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.caption = element_text(hjust = 0)
  )

ggsave(file.path(out_dir, "ifng_contrast_log2fc_scatter.png"), p_plain, width = 8, height = 8, dpi = 300)
ggsave(file.path(out_dir, "ifng_contrast_log2fc_scatter.pdf"), p_plain, width = 8, height = 8)
ggsave(file.path(out_dir, "ifng_contrast_log2fc_scatter_colored.png"), p_colored, width = 8, height = 8, dpi = 300)
ggsave(file.path(out_dir, "ifng_contrast_log2fc_scatter_colored.pdf"), p_colored, width = 8, height = 8)
