suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(ggplot2)
  library(readr)
})

in_deg <- "output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"
out_dir <- "output/DEA/plots/fgsea_jak_stat_interaction"
out_summary <- file.path(out_dir, "jak_stat_fgsea_summary.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

jak_stat_genes <- c(
  "APOL1", "APOL6", "IRF1", "CCL2", "CXCL9", "CXCL10",
  "ICAM1", "IDO1", "PDCD1LG2", "GBP1", "PSMB9", "CCDC68",
  "CTSS", "TMEM140", "TRIM21", "STAT1", "STAT2", "STAT3",
  "IRF2", "IRF8", "IRF9", "CXCL11"
)

pathways <- list(JAK_STAT_responsive = jak_stat_genes)

deg <- read_tsv(in_deg, show_col_types = FALSE)

make_ranks <- function(tbl, value_col) {
  value_sym <- rlang::sym(value_col)

  ranked_tbl <- tbl %>%
    filter(!is.na(genesymbol), genesymbol != "", !is.na(!!value_sym)) %>%
    group_by(genesymbol) %>%
    slice_max(order_by = abs(!!value_sym), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(!!value_sym), genesymbol)

  stats <- ranked_tbl[[value_col]]
  names(stats) <- ranked_tbl$genesymbol
  stats
}

run_fgsea_for_metric <- function(metric_col, metric_label) {
  ranks <- make_ranks(deg, metric_col)
  present_genes <- intersect(jak_stat_genes, names(ranks))

  fgsea_res <- fgsea(
    pathways = pathways,
    stats = ranks
  ) %>%
    as_tibble() %>%
    mutate(
      ranking_variable = metric_label,
      present_gene_count = length(present_genes),
      missing_genes = paste(setdiff(jak_stat_genes, present_genes), collapse = ","),
      leadingEdge = vapply(leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))
    ) %>%
    select(
      ranking_variable,
      pathway,
      size,
      ES,
      NES,
      pval,
      padj,
      leadingEdge,
      present_gene_count,
      missing_genes
    )

  plot_title <- paste("GSEA of JAK/STAT-responsive genes ranked by interaction", metric_label)
  plot_caption <- paste0(
    "Pathway genes present in ranking: ", length(present_genes), "/", length(jak_stat_genes)
  )

  p <- plotEnrichment(pathways[[1]], ranks) +
    labs(
      title = plot_title,
      x = metric_label,
      y = "Running enrichment score",
      caption = plot_caption
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0)
    )

  out_stub <- file.path(out_dir, paste0("jak_stat_fgsea_", metric_label))
  ggsave(paste0(out_stub, ".png"), p, width = 8, height = 5, dpi = 300)
  ggsave(paste0(out_stub, ".pdf"), p, width = 8, height = 5)

  fgsea_res
}

summary_tbl <- bind_rows(
  run_fgsea_for_metric("log2FoldChange", "log2fc"),
  run_fgsea_for_metric("stat", "stat")
)

write_tsv(summary_tbl, out_summary)
