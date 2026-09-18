suppressPackageStartupMessages({
  library(dplyr)
  library(fdrtool)
  library(readr)
  library(ggplot2)
})

if (exists("snakemake")) {
  in_deg <- snakemake@input[["in_deg"]]
  out_deg_cal <- snakemake@output[["out_deg_cal"]]
  out_summary <- snakemake@output[["out_summary"]]
  out_hist_raw_png <- snakemake@output[["out_hist_raw_png"]]
  out_hist_raw_pdf <- snakemake@output[["out_hist_raw_pdf"]]
  out_hist_cal_png <- snakemake@output[["out_hist_cal_png"]]
  out_hist_cal_pdf <- snakemake@output[["out_hist_cal_pdf"]]
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 7) {
    in_deg <- args[[1]]
    out_deg_cal <- args[[2]]
    out_summary <- args[[3]]
    out_hist_raw_png <- args[[4]]
    out_hist_raw_pdf <- args[[5]]
    out_hist_cal_png <- args[[6]]
    out_hist_cal_pdf <- args[[7]]
  } else {
    in_deg <- "output/DEA/objects/pyDESeq2_interaction/deg.tsv"
    out_deg_cal <- "output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"
    out_summary <- "output/DEA/objects/calibrate_pvalues_interaction/calibration_summary.tsv"
    out_hist_raw_png <- "output/DEA/plots/calibrate_pvalues_interaction/P_hist_raw.png"
    out_hist_raw_pdf <- "output/DEA/plots/calibrate_pvalues_interaction/P_hist_raw.pdf"
    out_hist_cal_png <- "output/DEA/plots/calibrate_pvalues_interaction/P_hist_calibrated.png"
    out_hist_cal_pdf <- "output/DEA/plots/calibrate_pvalues_interaction/P_hist_calibrated.pdf"
  }
}

dir.create(dirname(out_deg_cal), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_hist_raw_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_hist_cal_png), recursive = TRUE, showWarnings = FALSE)

deg <- read_tsv(in_deg, show_col_types = FALSE)
plot_subtitle <- if ("coefficient_tested" %in% names(deg)) {
  paste("Coefficient:", unique(deg$coefficient_tested)[1])
} else {
  "Coefficient: unknown"
}

calibrate <- function(res) {
  res_corrected <- res %>% filter(!is.na(stat), !is.na(pvalue))
  recalibrated <- fdrtool(res_corrected$stat, statistic = "normal", plot = FALSE, verbose = FALSE)
  res_corrected <- res_corrected %>%
    mutate(
      pvalue_raw = pvalue,
      padj_raw = padj,
      pvalue = recalibrated$pval,
      padj = p.adjust(recalibrated$pval, method = "BH"),
      null_sd = recalibrated$param[1, "sd"]
    )
  res_corrected
}

deg_cal <- calibrate(deg) %>% arrange(padj, pvalue)

write_tsv(deg_cal, out_deg_cal)

hist_raw <- ggplot(deg, aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), color = "black", fill = "grey70") +
  labs(title = "Raw P value histogram", subtitle = plot_subtitle)

hist_cal <- ggplot(deg_cal, aes(x = pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), color = "black", fill = "grey70") +
  labs(
    title = "Calibrated P value histogram",
    subtitle = paste0(plot_subtitle, ", null SD = ", round(unique(deg_cal$null_sd)[1], 3))
  )

ggsave(out_hist_raw_png, hist_raw, width = 6, height = 4, dpi = 300)
ggsave(out_hist_raw_pdf, hist_raw, width = 6, height = 4, dpi = 300)
ggsave(out_hist_cal_png, hist_cal, width = 6, height = 4, dpi = 300)
ggsave(out_hist_cal_pdf, hist_cal, width = 6, height = 4, dpi = 300)

summary_tab <- tibble(
  metric = c("raw_p_lt_0.05", "raw_padj_lt_0.05", "cal_p_lt_0.05", "cal_padj_lt_0.05", "null_sd"),
  value = c(
    sum(!is.na(deg$pvalue) & deg$pvalue < 0.05),
    sum(!is.na(deg$padj) & deg$padj < 0.05),
    sum(!is.na(deg_cal$pvalue) & deg_cal$pvalue < 0.05),
    sum(!is.na(deg_cal$padj) & deg_cal$padj < 0.05),
    unique(deg_cal$null_sd)[1]
  )
)

write_tsv(summary_tab, out_summary)
