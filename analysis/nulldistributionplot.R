#!/usr/bin/env Rscript
# Null distribution vs. observed statistic, for one pair, both directions
# (A-into-B and B-into-A).
#
# This reproduces the EXACT bootstrap procedure bitscoresumdifference's Go
# code uses (src/pipeline/reservoir.go's bootstrapMeans, src/stats/
# pairwise.go's empiricalPValue) -- not just a plot of raw reservoir
# scores. The null the p-value is actually tested against is a
# distribution of RESAMPLED MEDIANS (bootstrap draws of `resample_size`,
# median), which has much smaller spread than the raw individual scores
# in reservoir_matrix.tsv (basic CLT). Plotting raw scores next to the
# observed median would visually misrepresent how tight the actual test is
# -- so this script resamples in R exactly the way the Go code does,
# rather than shortcutting it.
#
# Requires: readr, dplyr, ggplot2
#   install.packages(c("readr", "dplyr", "ggplot2"))
#
# Inputs: real_matrix.tsv and reservoir_matrix.tsv as written by
# bitscoresumdifference's -dump-matrices flag (long format: sequence_id,
# origin, model, score).

library(readr)
library(dplyr)
library(ggplot2)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
real_matrix_path      <- "results/checkingbitsumdist/real_matrix.tsv" 
reservoir_matrix_path <- "results/checkingbitsumdist/reservoir_matrix.tsv"
model_a <- "SNORA17"
model_b <- "SNORA43"

# Should match (or be a deliberate multiple of) the actual run's
# -resample-size / -nshuffles if you want this plot to reproduce the
# exact p-value that run reported. resample_size defaults to -npool in
# the Go tool if unset there -- set it to whatever -npool was.
resample_size <- 100
n_boot        <- 200
seed          <- 42

output_path <- "null_vs_observed.svg"
stats_output_path <- "null_vs_observed_stats.txt"

# Manu "Kea" palette (github.com/G-Thomson/Manu) -- reused from earlier
# plots for visual consistency.
kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
                  "#D73202", "#272318", "#D3CDBF")

# --- load ----------------------------------------------------------------
real_df      <- read_tsv(real_matrix_path, show_col_types = FALSE)
reservoir_df <- read_tsv(reservoir_matrix_path, show_col_types = FALSE)

# --- helpers -------------------------------------------------------------

# median_score: median real score of `origin`'s pool against `target`.
median_score <- function(df, origin_name, target_name) {
  vals <- df %>% filter(origin == origin_name, model == target_name) %>% pull(score)
  if (length(vals) == 0) {
    stop(sprintf("no real_matrix rows for origin=%s, model=%s -- check model names", origin_name, target_name))
  }
  median(vals)
}

# bootstrap_null_stat: replicate reservoir.go's bootstrapMedians + the
# empiricalPValue construction exactly -- resample_size draws WITH
# REPLACEMENT from origin_name's reservoir scores against target_name,
# averaged, repeated n_boot times, each replicate's median then having the
# REAL self-baseline subtracted (fixed across all replicates, exactly as
# in the Go code) to form the null statistic distribution.
bootstrap_null_stat <- function(reservoir_df, origin_name, target_name, resample_size, n_boot, rng_seed) {
  pool_scores <- reservoir_df %>% filter(origin == origin_name, model == target_name) %>% pull(score)
  if (length(pool_scores) == 0) {
    stop(sprintf("no reservoir_matrix rows for origin=%s, model=%s -- check model names / that -dump-matrices was used",
                 origin_name, target_name))
  }
  if (resample_size > length(pool_scores)) {
    warning(sprintf("resample_size (%d) exceeds reservoir size for %s->%s (%d) -- results won't match the Go run's actual resampling",
                     resample_size, origin_name, target_name, length(pool_scores)))
  }
  set.seed(rng_seed)
  null_medians <- replicate(n_boot, median(sample(pool_scores, resample_size, replace = TRUE)))
  null_medians
}

# empirical_pvalue: same +1/+1 continuity-corrected one-sided test as
# Go's stats.empiricalPValue -- fraction of null replicates >= observed,
# bounded away from an overconfident p=0 by a finite bootstrap.
empirical_pvalue <- function(observed_stat, null_stat) {
  (sum(null_stat >= observed_stat) + 1) / (length(null_stat) + 1)
}

# --- compute both directions ----------------------------------------------
median_AA <- median_score(real_df, model_a, model_a)
median_AB <- median_score(real_df, model_a, model_b)
median_BB <- median_score(real_df, model_b, model_b)
median_BA <- median_score(real_df, model_b, model_a)

stat_a_into_b <- median_AB
stat_b_into_a <- median_BA

null_a_into_b <- bootstrap_null_stat(reservoir_df, model_a, model_b, resample_size, n_boot, seed)
null_b_into_a <- bootstrap_null_stat(reservoir_df, model_b, model_a, resample_size, n_boot, seed + 1)

p_a_into_b <- empirical_pvalue(stat_a_into_b, null_a_into_b)
p_b_into_a <- empirical_pvalue(stat_b_into_a, null_b_into_a)

plot_df <- bind_rows(
  data.frame(direction = sprintf("%s into %s", model_a, model_b), null_stat = null_a_into_b),
  data.frame(direction = sprintf("%s into %s", model_b, model_a), null_stat = null_b_into_a)
)
observed_df <- data.frame(
  direction = c(sprintf("%s into %s", model_a, model_b), sprintf("%s into %s", model_b, model_a)),
  observed  = c(stat_a_into_b, stat_b_into_a),
  pvalue    = c(p_a_into_b, p_b_into_a)
)
observed_df$label <- sprintf("observed = %.2f, p = %.4g", observed_df$observed, observed_df$pvalue)

# --- plot ------------------------------------------------------------
p <- ggplot(plot_df, aes(x = null_stat)) +
  geom_density(fill = kea_palette[1], alpha = 0.85, color = "white", linewidth = 0.2) +
  geom_vline(data = observed_df, aes(xintercept = observed), color = kea_palette[5], linewidth = 0.9) +
  geom_text(data = observed_df, aes(x = observed, y = Inf, label = label),
            color = kea_palette[5], angle = 90, vjust = -0.6, hjust = 1.05, size = 3.3) +
  facet_wrap(~ direction, ncol = 1, scales = "free") +
  labs(
    x = expression(paste("null statistic (bootstrap resampled ", median(shuffled), ")")),
    y = "count"
  ) +
  theme_publication(base_size = 12, legend_position = "right") +
  theme(strip.text = element_text(face = "bold"))

#print(p)
ggsave(output_path, p, width = 8.5, height = 6.5, dpi = 300, bg = "white")

stats_summary <- list(
  model_a = model_a,
  model_b = model_b,
  resample_size = resample_size,
  n_boot = n_boot,
  observed_df = observed_df
)

lines <- c(
  sprintf("model_a=%s", stats_summary$model_a),
  sprintf("model_b=%s", stats_summary$model_b),
  sprintf("resample_size=%d", stats_summary$resample_size),
  sprintf("n_boot=%d", stats_summary$n_boot)
)
for (i in seq_len(nrow(stats_summary$observed_df))) {
  row <- stats_summary$observed_df[i, ]
  lines <- c(
    lines,
    sprintf("direction=%s", row$direction),
    sprintf("observed=%.6g", row$observed),
    sprintf("p_value=%.6g", row$pvalue)
  )
}
writeLines(lines, con = stats_output_path)

cat("\n")
print(observed_df[, c("direction", "observed", "pvalue")])