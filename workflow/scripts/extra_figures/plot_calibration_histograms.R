suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

if (exists("snakemake")) {
  in_raw <- snakemake@input[["raw"]]
  in_calibrated <- snakemake@input[["calibrated"]]
  out_csv <- snakemake@output[["csv"]]
  out_png <- snakemake@output[["png"]]
  out_pdf <- snakemake@output[["pdf"]]
} else {
  in_raw <- "results/dea/interaction.tsv"
  in_calibrated <- "results/dea/interaction_calibrated.tsv"
  out_csv <- "results/extra_figures/source_data/calibration_histograms.tsv"
  out_png <- "plots/extra_figures/calibration_histograms.png"
  out_pdf <- "plots/extra_figures/calibration_histograms.pdf"
}

for (path in c(out_csv, out_png, out_pdf)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

raw <- read_tsv(in_raw, show_col_types = FALSE) %>%
  select(gene_id, raw = pvalue)
calibrated <- read_tsv(in_calibrated, show_col_types = FALSE) %>%
  select(gene_id, calibrated = pvalue)
source_data <- inner_join(raw, calibrated, by = "gene_id") %>%
  pivot_longer(c(raw, calibrated), names_to = "distribution", values_to = "pvalue") %>%
  filter(!is.na(pvalue))
write_tsv(source_data, out_csv)

plot <- ggplot(source_data, aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), boundary = 0,
                 fill = "#4C78A8", colour = "white") +
  facet_wrap(~distribution, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_classic(base_family = "Arial", base_size = 8) +
  labs(x = "P-value", y = "Genes", title = "Empirical-null calibration")

ggsave(out_png, plot, width = 7, height = 7, dpi = 300)
ggsave(out_pdf, plot, width = 7, height = 7, device = cairo_pdf)
