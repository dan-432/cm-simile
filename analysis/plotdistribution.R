#!/usr/bin/env Rscript
# Probability distribution of sequence scores (from cmsearch) between two models (cross and against themselves)
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
library(ggplot2)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
input_path  <- "results/checkingbitsumdist/real_matrix.tsv"             # bitscoresumdifference output
null_path <- "results/checkingbitsumdist/reservoir_matrix.tsv"             # bitscoresumdifference output

modelA <- "SNORA17"
modelB <- "SNORA43" # "Pyrobac-1" #

output_path <- sprintf("%s_%s_distribution.pdf", modelA, modelB)
stats_output_path <- sprintf("%s_%s_distribution_stats.txt", modelA, modelB)

# Manu "Kea" palette (github.com/G-Thomson/Manu)
kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
                  "#D73202", "#272318", "#D3CDBF")

# has columns sequence_id	origin	model	score
scores <- read_tsv(input_path, show_col_types = FALSE)
null_scores <- read_tsv(null_path, show_col_types = FALSE)

# limit number of null scores to same size as the real scores for plotting purposes (otherwise the null distribution is too dense to see the real distribution)
null_scores <- null_scores[sample(nrow(null_scores), min(nrow(scores), nrow(null_scores))), ]

# Convert tibble column slices to plain numeric vectors so ggplot/ks.test
# receive numeric vectors rather than 1-column tibbles (which are lists).
self_scores_a <- scores$score[scores$model == modelA & scores$origin == modelA]
self_scores_b <- scores$score[scores$model == modelB & scores$origin == modelB]
b_scored_on_model_a <- scores$score[scores$model == modelA & scores$origin == modelB]
a_scored_on_model_b <- scores$score[scores$model == modelB & scores$origin == modelA]
a_null_scores <- null_scores$score[null_scores$model == modelA & null_scores$origin == modelB]
b_null_scores <- null_scores$score[null_scores$model == modelB & null_scores$origin == modelA]

## --- plot ------------------------------------------------------------
p <- ggplot() +
  geom_density(aes(x = self_scores_a, color = "self_scores_a"), linewidth = 3, alpha = 0.7) +
  geom_density(aes(x = self_scores_b, color = "self_scores_b"), linewidth = 3, alpha = 0.7) +
  geom_density(aes(x = b_scored_on_model_a, color = "b_scored_on_model_a"), linewidth = 3) +
  geom_density(aes(x = a_scored_on_model_b, color = "a_scored_on_model_b"), linewidth = 3) +
  geom_density(aes(x = a_null_scores, color = "a_null_scores"), linewidth = 3) +
  geom_density(aes(x = b_null_scores, color = "b_null_scores"), linewidth = 3) +
  scale_color_manual(values = c(
    self_scores_a = kea_palette[1],
    self_scores_b = kea_palette[3],
    b_scored_on_model_a = kea_palette[2],
    a_scored_on_model_b = kea_palette[4],
    a_null_scores = kea_palette[6],
    b_null_scores = kea_palette[7]
  ), labels = c(
    self_scores_a = sprintf("%s self scores", modelA),
    self_scores_b = sprintf("%s self scores", modelB),
    b_scored_on_model_a = sprintf("%s scored on %s", modelB, modelA),
    a_scored_on_model_b = sprintf("%s scored on %s", modelA, modelB),
    a_null_scores = sprintf("%s null scores", modelA),
    b_null_scores = sprintf("%s null scores", modelB)
  )) +
  labs(
    x = "Bit Score",
    y = "Density",
    color = NULL
  ) +
  theme_publication(base_size = 36, legend_position = "bottom") +
  theme(
    axis.line = element_line(linewidth = 1.8),
    axis.ticks = element_line(linewidth = 1.5),
    legend.key.size = unit(1.65, "cm"),
    legend.key = element_blank(),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(
    nrow = 3,
    byrow = TRUE,
    override.aes = list(linewidth = 3, alpha = 1)
  ))

# Save the plot
ggsave(output_path, p, width =16, height = 16, dpi = 300, bg = "white")


# do a KS test between cross scores and null scores
ks_test_b_on_a <- ks.test(b_scored_on_model_a, a_null_scores)
ks_test_a_on_b <- ks.test(a_scored_on_model_b, b_null_scores)

# Print the KS test results
cat(sprintf("KS test for %s scored on %s vs null: D = %.4f, p-value = %.4g\n", modelB, modelA, ks_test_b_on_a$statistic, ks_test_b_on_a$p.value))
cat(sprintf("KS test for %s scored on %s vs null: D = %.4f, p-value = %.4g\n", modelA, modelB, ks_test_a_on_b$statistic, ks_test_a_on_b$p.value))

# do ks test between cross scores and self scores
ks_test_b_on_a_self <- ks.test(b_scored_on_model_a, self_scores_a)
ks_test_a_on_b_self <- ks.test(a_scored_on_model_b, self_scores_b)

# Print the KS test results for cross vs self scores
cat(sprintf("KS test for %s scored on %s vs self: D = %.4f, p-value = %.4g\n", modelB, modelA, ks_test_b_on_a_self$statistic, ks_test_b_on_a_self$p.value))
cat(sprintf("KS test for %s scored on %s vs self: D = %.4f, p-value = %.4g\n", modelA, modelB, ks_test_a_on_b_self$statistic, ks_test_a_on_b_self$p.value))

stats_summary <- list(
  model_a = modelA,
  model_b = modelB,
  ks_b_on_a_vs_null_d = ks_test_b_on_a$statistic,
  ks_b_on_a_vs_null_p = ks_test_b_on_a$p.value,
  ks_a_on_b_vs_null_d = ks_test_a_on_b$statistic,
  ks_a_on_b_vs_null_p = ks_test_a_on_b$p.value,
  ks_b_on_a_vs_self_d = ks_test_b_on_a_self$statistic,
  ks_b_on_a_vs_self_p = ks_test_b_on_a_self$p.value,
  ks_a_on_b_vs_self_d = ks_test_a_on_b_self$statistic,
  ks_a_on_b_vs_self_p = ks_test_a_on_b_self$p.value
)

lines <- c(
  sprintf("model_a=%s", stats_summary$model_a),
  sprintf("model_b=%s", stats_summary$model_b),
  sprintf("ks_b_on_a_vs_null_d=%.6g", stats_summary$ks_b_on_a_vs_null_d),
  sprintf("ks_b_on_a_vs_null_p=%.6g", stats_summary$ks_b_on_a_vs_null_p),
  sprintf("ks_a_on_b_vs_null_d=%.6g", stats_summary$ks_a_on_b_vs_null_d),
  sprintf("ks_a_on_b_vs_null_p=%.6g", stats_summary$ks_a_on_b_vs_null_p),
  sprintf("ks_b_on_a_vs_self_d=%.6g", stats_summary$ks_b_on_a_vs_self_d),
  sprintf("ks_b_on_a_vs_self_p=%.6g", stats_summary$ks_b_on_a_vs_self_p),
  sprintf("ks_a_on_b_vs_self_d=%.6g", stats_summary$ks_a_on_b_vs_self_d),
  sprintf("ks_a_on_b_vs_self_p=%.6g", stats_summary$ks_a_on_b_vs_self_p)
)
writeLines(lines, con = stats_output_path)

# plot culmulative distribution functions (CDFs) for null, cross scores and self scores
cdf_self_plot <- ggplot() +
  stat_ecdf(aes(x = b_scored_on_model_a, color = "b_scored_on_model_a"), linewidth = 3) +
  stat_ecdf(aes(x = a_scored_on_model_b, color = "a_scored_on_model_b"), linewidth = 3) +
  stat_ecdf(aes(x = self_scores_a, color = "self_scores_a"), linewidth = 3) +
  stat_ecdf(aes(x = self_scores_b, color = "self_scores_b"), linewidth = 3) +
  stat_ecdf(aes(x = a_null_scores, color = "a_null_scores"), linewidth = 3) +
  stat_ecdf(aes(x = b_null_scores, color = "b_null_scores"), linewidth = 3) +
  scale_color_manual(values = c(
    b_scored_on_model_a = kea_palette[2],
    a_scored_on_model_b = kea_palette[4],
    self_scores_a = kea_palette[1],
    self_scores_b = kea_palette[3],
    a_null_scores = kea_palette[6],
    b_null_scores = kea_palette[7]
  ), labels = c(
    b_scored_on_model_a = sprintf("%s scored on %s", modelB, modelA),
    a_scored_on_model_b = sprintf("%s scored on %s", modelA, modelB),
    self_scores_a = sprintf("%s self scores", modelA),
    self_scores_b = sprintf("%s self scores", modelB),
    a_null_scores = sprintf("%s null scores", modelA),
    b_null_scores = sprintf("%s null scores", modelB)
  )) +
  labs(
    x = "Bit Score",
    y = "Cumulative Distribution Function (CDF)",
    color = NULL
  ) +
  theme_publication(base_size = 36, legend_position = "bottom") +
  theme(
    axis.line = element_line(linewidth = 1.8),
    axis.ticks = element_line(linewidth = 1.5),
    legend.key.size = unit(1.65, "cm"),
    legend.key = element_blank(),
    legend.key.width = unit(2.5, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.box = "horizontal"
  ) +
  guides(color = guide_legend(
    nrow = 3,
    byrow = TRUE,
    override.aes = list(linewidth = 3, alpha = 1)
  ))

# Save the CDF plot for cross vs self scores
cdf_self_output_path <- sprintf("%s_%s_cdf_cross_vs_self.pdf", modelA, modelB)
ggsave(cdf_self_output_path, cdf_self_plot, width = 16, height = 16, dpi = 300, bg = "white")