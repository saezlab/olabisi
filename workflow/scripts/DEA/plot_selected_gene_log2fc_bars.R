suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

case_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_case__IFNGpos_vs_IFNGneg/deg.tsv"
control_path <- "output/DEA/objects/pyDESeq2_byGrp/HR_control__IFNGpos_vs_IFNGneg/deg.tsv"
out_dir <- "output/DEA/plots/log2fc_barplot_selected_genes"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selected_genes <- c(
  "APOL1", "APOL6", "IRF1", "CCL2", "CXCL9", "CXCL10",
  "ICAM1", "IDO1", "PDCD1LG2", "GBP1", "PSMB9", "CCDC68",
  "CTSS", "TMEM140", "TRIM21", "STAT1", "STAT2", "STAT3",
  "IRF2", "IRF8", "IRF9", "CXCL11"
)

read_deg <- function(path, contrast_label) {
  read_tsv(path, show_col_types = FALSE) %>%
    select(genesymbol, log2FoldChange, padj) %>%
    mutate(contrast = contrast_label)
}

case_deg <- read_deg(case_path, "HR case")
control_deg <- read_deg(control_path, "HR control")

plot_df <- bind_rows(case_deg, control_deg) %>%
  filter(genesymbol %in% selected_genes) %>%
  mutate(
    genesymbol = factor(genesymbol, levels = selected_genes),
    contrast = factor(contrast, levels = c("HR case", "HR control")),
    sig_label = if_else(!is.na(padj) & padj < 0.05, "*", ""),
    sign_dir = if_else(log2FoldChange >= 0, 1, -1)
  )

missing_pairs <- expand.grid(
  genesymbol = selected_genes,
  contrast = c("HR case", "HR control"),
  stringsAsFactors = FALSE
) %>%
  anti_join(
    plot_df %>% transmute(genesymbol = as.character(genesymbol), contrast = as.character(contrast)),
    by = c("genesymbol", "contrast")
  )

if (nrow(missing_pairs) > 0) {
  warning(
    "Missing selected genes in one or more contrasts: ",
    paste(paste0(missing_pairs$genesymbol, " [", missing_pairs$contrast, "]"), collapse = ", ")
  )
}

facet_limits <- plot_df %>%
  group_by(genesymbol) %>%
  summarise(
    ymin = min(log2FoldChange, 0, na.rm = TRUE),
    ymax = max(log2FoldChange, 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    span = pmax(ymax - ymin, 0.5),
    star_y = if_else(ymax >= 0, ymax + 0.10 * span, ymin - 0.10 * span)
  )

plot_df <- plot_df %>%
  left_join(facet_limits %>% select(genesymbol, star_y), by = "genesymbol") %>%
  mutate(star_vjust = if_else(log2FoldChange >= 0, 0, 1))

p <- ggplot(plot_df, aes(x = contrast, y = log2FoldChange, fill = contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey45") +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  geom_text(
    data = filter(plot_df, sig_label != ""),
    aes(y = star_y, label = sig_label),
    size = 6,
    vjust = filter(plot_df, sig_label != "")$star_vjust
  ) +
  facet_wrap(~ genesymbol, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("HR case" = "#b2182b", "HR control" = "#313695")) +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title = "log2FC for selected genes in IFNG+ vs IFNG- contrasts by disease group",
    x = NULL,
    y = "log2FC",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "top"
  )

ggsave(file.path(out_dir, "selected_genes_log2fc_barplot.png"), p, width = 14, height = 10, dpi = 300)
ggsave(file.path(out_dir, "selected_genes_log2fc_barplot.pdf"), p, width = 14, height = 10)
