suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(ggplot2)
  library(readr)
})

if (exists("snakemake")) {
  in_degs <- snakemake@input[["degs"]]
  in_network <- snakemake@input[["network"]]
  out_csv <- snakemake@output[["csv"]]
  out_png <- snakemake@output[["png"]]
  out_pdf <- snakemake@output[["pdf"]]
  seed <- as.integer(snakemake@params[["seed"]])
} else {
  in_degs <- "results/dea/interaction_calibrated.tsv"
  in_network <- "data/prior_knowledge/progeny_20260610.tsv"
  out_csv <- "results/extra_figures/source_data/jak_stat_fgsea.tsv"
  out_png <- "plots/extra_figures/jak_stat_fgsea.png"
  out_pdf <- "plots/extra_figures/jak_stat_fgsea.pdf"
  seed <- 0L
}

for (path in c(out_csv, out_png, out_pdf)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
set.seed(seed)

degs <- read_tsv(in_degs, show_col_types = FALSE) %>%
  filter(!is.na(genesymbol), !is.na(stat)) %>%
  group_by(genesymbol) %>%
  filter(n() == 1L) %>%
  ungroup()
ranks <- degs$stat
names(ranks) <- degs$genesymbol
ranks <- sort(ranks, decreasing = TRUE)

network <- read_tsv(in_network, show_col_types = FALSE)
genes <- network %>%
  filter(source == "JAK-STAT", target %in% names(ranks)) %>%
  pull(target) %>%
  unique()
if (length(genes) < 10L) stop("Too few JAK-STAT targets overlap the ranked DE table")

result <- fgsea(
  pathways = list(`JAK-STAT` = genes),
  stats = ranks,
  minSize = 10,
  maxSize = Inf,
  eps = 0
) %>%
  as.data.frame() %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))
write_tsv(result, out_csv)

plot <- plotEnrichment(genes, ranks) +
  theme_classic(base_family = "Arial", base_size = 8) +
  labs(title = "JAK-STAT interaction-statistic FGSEA")
ggsave(out_png, plot, width = 7, height = 5, dpi = 300)
ggsave(out_pdf, plot, width = 7, height = 5, device = cairo_pdf)
