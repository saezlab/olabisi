if (exists("snakemake")) {
  out_summary <- snakemake@output[["summary"]]
} else {
  out_summary <- "results/provenance/renv_packages.txt"
}

dir.create(dirname(out_summary), recursive = TRUE, showWarnings = FALSE)
packages <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
packages <- packages[order(packages$Package), , drop = FALSE]
writeLines(c(capture.output(sessionInfo()), "", capture.output(packages)), out_summary)
