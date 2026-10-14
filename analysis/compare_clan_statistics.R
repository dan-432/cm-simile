#!/usr/bin/env Rscript
# Compare the three pairwise similarity approaches by how well they separate
# within-clan pairs from pairs in different clans.
#
# The benchmark populations are analyzed independently because the three runs
# were made with different population directories. The reported ANOVA values
# are useful for ranking the statistics, but pairwise family comparisons are
# not fully independent; use the effect sizes and AUC alongside the p-values.
#
# Requires: readr
#   install.packages("readr")

library(readr)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
benchmark_inputs <- list(
  bitscoresumdifference = list(
    path = "results/bitscoresum_benchmark_5687123/raw_results.tsv",
    clans = "results/clan_population_5687123/population_clans.tsv"
  ),
  ksmwcompare = list(
    path = "results/ksmwcompare_benchmark_5740367/raw_results.tsv",
    clans = "results/clan_population_5740367/population_clans.tsv"
  ),
  minmaxcmcompare = list(
    path = "results/minmax_benchmark_5664258/raw_results.tsv",
    clans = "results/clan_population_5664258/population_clans.tsv"
  )
)

output_path <- file.path("analysis", "clan_statistic_anova.tsv")

# Convert each method's native output into a common score where larger means
# stronger evidence of similarity, making effect sizes and AUC comparable.
method_definitions <- list(
  bitscoresumdifference = list(
    columns = c("stat_a_into_b", "stat_b_into_a"),
    score = function(df) -(df$stat_a_into_b + df$stat_b_into_a)
  ),
  ksmwcompare = list(
    columns = c("ks_d_a_into_b", "ks_d_b_into_a"),
    score = function(df) 1 - (df$ks_d_a_into_b + df$ks_d_b_into_a) / 2
  ),
  minmaxcmcompare = list(
    columns = "similarity_ratio",
    score = function(df) df$similarity_ratio
  )
)

check_columns <- function(df, required, path) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing column(s) %s in %s. Available columns: %s",
      paste(missing, collapse = ", "), path, paste(names(df), collapse = ", ")
    ))
  }
}

calculate_auc <- function(score, within_clan) {
  positive <- within_clan == "within_clan"
  n_positive <- sum(positive)
  n_negative <- sum(!positive)
  if (n_positive == 0 || n_negative == 0) {
    return(NA_real_)
  }

  score_ranks <- rank(score, ties.method = "average")
  (sum(score_ranks[positive]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

analyze_method <- function(method_name, config, definition) {
  df <- read_tsv(config$path, show_col_types = FALSE)
  clans <- read_tsv(config$clans, show_col_types = FALSE)
  check_columns(df, c("model_a", "model_b", definition$columns), config$path)

  df <- prepare_clan_comparison(
    df,
    stat_col = definition$columns[1],
    clans = clans,
    accession_col = "accession",
    clan_col = "clan",
    model_a_col = "model_a",
    model_b_col = "model_b"
  )
  df$similarity <- definition$score(df)
  df <- df[is.finite(df$similarity) & !is.na(df$comparison), ]

  if (nrow(df) < 3 || any(table(df$comparison) < 2)) {
    stop(sprintf("Not enough observations in %s for both comparison groups.", method_name))
  }

  fit <- aov(similarity ~ comparison, data = df)
  anova_table <- summary(fit)[[1]]
  between_ss <- anova_table["comparison", "Sum Sq"]
  total_ss <- sum(anova_table[, "Sum Sq"])
  within_values <- df$similarity[df$comparison == "within_clan"]
  other_values <- df$similarity[df$comparison == "everything_else"]
  pooled_sd <- sqrt(((length(within_values) - 1) * stats::var(within_values) +
    (length(other_values) - 1) * stats::var(other_values)) /
    (length(within_values) + length(other_values) - 2))

  data.frame(
    method = method_name,
    n_total = nrow(df),
    n_within_clan = length(within_values),
    n_everything_else = length(other_values),
    mean_within_clan = mean(within_values),
    mean_everything_else = mean(other_values),
    median_within_clan = median(within_values),
    median_everything_else = median(other_values),
    anova_f = unname(anova_table["comparison", "F value"]),
    anova_p = unname(anova_table["comparison", "Pr(>F)"]),
    eta_squared = between_ss / total_ss,
    standardized_mean_difference = (mean(within_values) - mean(other_values)) / pooled_sd,
    within_clan_auc = calculate_auc(df$similarity, df$comparison),
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  lapply(names(benchmark_inputs), function(method_name) {
    analyze_method(method_name, benchmark_inputs[[method_name]], method_definitions[[method_name]])
  })
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(results, output_path)
print(results)
cat(sprintf("Wrote %s\n", output_path))
