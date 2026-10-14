#!/usr/bin/env Rscript
# Heatmap of pairwise family comparisons
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
input_path <- "results/bitscoresum_benchmark_5687123/raw_results.tsv"
input_clan_path <- "results/clan_population_5687123/population_clans.tsv"
output_path <- file.path("analysis", "bitscoresum_pvalue_heatmap.pdf")

df <- read_tsv(input_path, show_col_types = FALSE)
clans <- read_tsv(input_clan_path, show_col_types = FALSE)

required_pvalue_cols <- c("pvalue_a_into_b", "pvalue_b_into_a")
missing_pvalue_cols <- setdiff(required_pvalue_cols, names(df))
if (length(missing_pvalue_cols) > 0) {
  stop(sprintf(
    "Missing required p-value column(s): %s",
    paste(missing_pvalue_cols, collapse = ", ")
  ))
}

# Larger values represent stronger evidence, so plot p-values as -log10(p).
# Clamp at the smallest representable p-value to avoid infinite values.
min_pvalue <- .Machine$double.xmin
df$minus_log10_p_a_into_b <- -log10(pmax(df$pvalue_a_into_b, min_pvalue))
df$minus_log10_p_b_into_a <- -log10(pmax(df$pvalue_b_into_a, min_pvalue))

df <- prepare_clan_comparison(
  df,
  stat_col = "minus_log10_p_a_into_b",
  clans = clans,
  accession_col = "accession",
  clan_col = "clan",
  model_a_col = "model_a",
  model_b_col = "model_b"
)

# --- plot ------------------------------------------------------------
result <- plot_heatmap(
  df,
  group_col = "comparison",
  stat_col = "minus_log10_p_a_into_b",
  reverse_stat_col = "minus_log10_p_b_into_a",
  order_group_col = "clan",
  model_a_col = "model_a",
  model_b_col = "model_b",
  title = expression(-log[10](p)),
  save_path = output_path,
  width = 9,
  height = 9,
  dpi = 300
)