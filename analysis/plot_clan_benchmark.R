#!/usr/bin/env Rscript
# Violin-and-point plot of the combined bidirectional KS statistic for
# within-clan pairs versus other pairs.
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
library(ggplot2)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
input_path <- "results/ksmwcompare_benchmark_5740367/raw_results.tsv"
input_clan_path <- "results/clan_population_5740367/population_clans.tsv"
output_path <- file.path("analysis", "ks_ratio_violinplot.pdf")
stats_output_path <- file.path("analysis", "ks_ratio_violinplot_stats.txt")

stat_column <- "ks_d_sum"
stat_components <- c("ks_d_a_into_b", "ks_d_b_into_a")

# Manu "Kea" palette (github.com/G-Thomson/Manu)
kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
                  "#D73202", "#272318", "#D3CDBF")

# --- load & label --------------------------------------------------------
df <- read_tsv(input_path, show_col_types = FALSE)
clans <- read_tsv(input_clan_path, show_col_types = FALSE)

required_cols <- c("model_a", "model_b", stat_components)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "missing column(s) %s in %s -- available columns: %s",
    paste(missing_cols, collapse = ", "), input_path, paste(names(df), collapse = ", ")
  ))
}
df[[stat_column]] <- df[[stat_components[1]]] + df[[stat_components[2]]]
top_n <- max(1L, ceiling(nrow(df) * 0.03))
top_idx <- order(df[[stat_column]], decreasing = TRUE)[seq_len(top_n)]
df$pair_label <- NA_character_
df$pair_label[top_idx] <- paste(df$model_a[top_idx], df$model_b[top_idx], sep = " / ")

df <- prepare_clan_comparison(
  df,
  stat_col = stat_column,
  clans = clans,
  accession_col = "accession",
  clan_col = "clan",
  model_a_col = "model_a",
  model_b_col = "model_b"
)

# --- plot ------------------------------------------------------------
plot_result <- plot_clan_box_dotplot(
  df,
  stat_col = stat_column,
  y_label = "KS statistic (bidirectional sum)",
  title = NULL,
  save_path = output_path,
  stats_output_path = stats_output_path,
  width = 6,
  height = 5,
  fill_palette = c(within_clan = kea_palette[1], everything_else = kea_palette[2]),
  label_col = "pair_label"
)

print(plot_result$plot)


# ---- modeling ----

# df$within_clan <- as.integer(
#   !is.na(df$clan_a) &
#   !is.na(df$clan_b) &
#   df$clan_a != "" &
#   df$clan_a == df$clan_b
# )

# library(lme4)

# m1 <- lmer(
#   median_score_a_on_b ~ within_clan +
#     (1 | model_a) +
#     (1 | model_b),
#   data = df
# )

# ##summary(m1)
# fixef(m1)

# df$cross_sum <-
#   df$median_score_a_on_b +
#   df$median_score_b_on_a

# m2 <- lmer(
#   cross_sum ~ within_clan +
#     (1 | model_a) +
#     (1 | model_b),
#   data = df
# )

# ##summary(m2)
# fixef(m2)

# m3 <- lmer(
#   stat_a_into_b ~ within_clan +
#     (1 | model_a) +
#     (1 | model_b),
#   data = df
# )

# #summary(m3)
# fixef(m3)

# library(ggplot2)

# ggplot(df,
#        aes(median_score_b_on_b,
#            median_score_a_on_b,
#            color = factor(within_clan))) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm")


# df$cross_sum <-
#     df$median_score_a_on_b +
#     df$median_score_b_on_a

# df$self_sum <-
#     df$median_score_a_on_a +
#     df$median_score_b_on_b

# m4 <- lm(
#   cross_sum ~ self_sum + within_clan,
#   data = df
# )
# summary(m4)

# m_cross <- lm(
#     cross_sum ~ within_clan,
#     data=df
# )

# summary(m_cross)

# anova(m_cross, m4)



# # library(readr)
# # library(ggplot2)

# # # --- config ------------------------------------------------------------
# # input_path  <- "results/minmax_benchmark_5664258/raw_results.tsv"             # bitscoresumdifference/clan-benchmark output
# # input_clan_path  <- "results/clan_population_5664258/population_clans.tsv"
# # output_path <- "minmax_boxplot.png"

# # # Manu "Kea" palette (github.com/G-Thomson/Manu)
# # kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
# #                   "#D73202", "#272318", "#D3CDBF")

# # # --- load & label --------------------------------------------------------
# # df <- read_tsv(input_path, show_col_types = FALSE)
# # clans <- read_tsv(input_clan_path, show_col_types = FALSE)

# # # Build a lookup of accession -> clan, then annotate each pair.
# # clan_map <- setNames(clans$clan, clans$accession)
# # df$clan_a <- clan_map[df$model_a]
# # df$clan_b <- clan_map[df$model_b]

# # df$comparison <- ifelse(
# #   !is.na(df$clan_a) & !is.na(df$clan_b) & df$clan_a != "" & df$clan_a == df$clan_b,
# #   "within_clan",
# #   "everything_else"
# # )

# # df$comparison <- factor(df$comparison, levels = c("within_clan", "everything_else"))

# # # quick Wilcoxon (Mann-Whitney) test
# # wt <- wilcox.test(similarity_ratio ~ comparison, data = df)
# # n_labels <- table(df$comparison)
# # subtitle_txt <- sprintf(
# #   "within_clan n=%d, everything_else n=%d  |  Wilcoxon p = %.3g",
# #   n_labels["within_clan"], n_labels["everything_else"], wt$p.value
# # )

# # # --- plot ------------------------------------------------------------
# # p <- ggplot(df, aes(x = comparison, y = similarity_ratio, fill = comparison)) +
# #   geom_boxplot(outlier.shape = 21, outlier.alpha = 0.6, width = 0.5) +
# #   scale_fill_manual(values = c(within_clan = kea_palette[1],
# #                                 everything_else = kea_palette[2])) +
# #   labs(
# #     x = NULL,
# #     y = "Similarity Ratio",
# #     title = "Cross-recognition statistic: within-clan vs. everything else",
# #     subtitle = subtitle_txt,
# #     fill = NULL
# #   ) +
# #   theme_minimal(base_size = 13) +
# #   theme(legend.position = "none")

# # print(p)
# # ggsave(output_path, p, width = 6, height = 5, dpi = 300)