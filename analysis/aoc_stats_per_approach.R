#!/usr/bin/env Rscript
# AUC and Wilcoxon summary for each pairwise scoring approach

library(readr)
source(file.path("analysis", "publication_plot_helpers.R"))

clans <- read_tsv("results/clan_population_5687123/population_clans.tsv",
                  show_col_types = FALSE)

# One entry per approach: results file, and how to build its per-pair statistic.
# Fill in the paths and column names to match your raw_results files.
approaches <- list(
  list(name = "Common ancestor (similarity ratio)",
       path = "results/minmax_benchmark_5664258/raw_results.tsv",
       stat = function(d) d$similarity_ratio),
  list(name = "Median bit score vs. null",
       path = "results/bitscoresum_benchmark_5687123/raw_results.tsv",
       stat = function(d) d$median_score_a_on_b + d$median_score_b_on_a),
  list(name = "KS distribution vs. null",
       path = "results/ksmwcompare_benchmark_5740367/raw_results.tsv",
       stat = function(d) d$ks_d_a_into_b + d$ks_d_b_into_a)
)

for (a in approaches) {
  df <- read_tsv(a$path, show_col_types = FALSE)
  df$stat <- a$stat(df)
  df <- prepare_clan_comparison(df, stat_col = "stat", clans = clans)
  s <- make_comparison_summary(df, stat_col = "stat")
  cat(sprintf("%s\n  %s\n  %s\n  %s\n\n",
              a$name, s$n_text, s$auc_text, s$p_text))
}