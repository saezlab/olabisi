suppressPackageStartupMessages({
  library(dplyr)
  library(fdrtool)
  library(readr)
})

if (exists("snakemake")) {
  in_degs <- snakemake@input[["degs"]]
  out_degs <- snakemake@output[["degs"]]
  out_summary <- snakemake@output[["summary"]]
} else {
  in_degs <- "results/dea/interaction.tsv"
  out_degs <- "results/dea/interaction_calibrated.tsv"
  out_summary <- "results/dea/calibration_summary.tsv"
}

dir.create(dirname(out_degs), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)

degs <- read_tsv(in_degs, show_col_types = FALSE)
valid <- degs %>% filter(!is.na(stat), !is.na(pvalue))
fit <- fdrtool(valid$stat, statistic = "normal", plot = FALSE, verbose = FALSE)
calibrated <- valid %>%
  mutate(
    pvalue_raw = pvalue,
    padj_raw = padj,
    pvalue = fit$pval,
    padj = p.adjust(fit$pval, method = "BH"),
    null_sd = fit$param[1, "sd"]
  ) %>%
  arrange(padj, pvalue, gene_id)

summary_tab <- tibble(
  metric = c(
    "rows_input", "rows_calibrated", "raw_p_lt_0.05", "raw_padj_lt_0.05",
    "cal_p_lt_0.05", "cal_padj_lt_0.05", "null_sd"
  ),
  value = c(
    nrow(degs), nrow(calibrated),
    sum(!is.na(degs$pvalue) & degs$pvalue < 0.05),
    sum(!is.na(degs$padj) & degs$padj < 0.05),
    sum(calibrated$pvalue < 0.05),
    sum(calibrated$padj < 0.05),
    unique(calibrated$null_sd)[1]
  )
)

write_tsv(calibrated, out_degs)
write_tsv(summary_tab, out_summary)
