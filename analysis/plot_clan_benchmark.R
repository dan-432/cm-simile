#!/usr/bin/env Rscript
# Boxplot of stat_a_into_b (now median_b_on_a per Dan's redefinition),
# within-clan pairs vs. everything else -- bitscoresumdifference's
# analogue of the minmax_boxplot.R script this was based on.
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

#!/usr/bin/env Rscript
# Boxplot of stat_a_into_b + stat_b_into_a (a combined, symmetric
# statistic), within-clan pairs vs. everything else --
# bitscoresumdifference's analogue of the minmax_boxplot.R script this
# was based on. Both stat_* columns now hold medians rather than means
# (column names unchanged from the original design), and this script
# sums both directions rather than plotting one alone.
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
library(ggplot2)

# --- config ------------------------------------------------------------
input_path      <- "results/bitscoresum_benchmark_5687123/raw_results.tsv"  # bitscoresumdifference output
input_clan_path <- "results/clan_population_5687123/population_clans.tsv"
output_path     <- "bitscoresumdifference_boxplot.svg"

# Combined symmetric statistic: stat_a_into_b + stat_b_into_a. Column
# NAMES are unchanged from the original mean-based design (still
# "stat_a_into_b"/"stat_b_into_a", "mean_score_*" etc.) even though the
# underlying values are now medians, not means -- Dan confirmed this
# directly, so no renamed-column guessing needed here, unlike the
# previous version of this script.
stat_column <- "combined_stat"

# Manu "Kea" palette (github.com/G-Thomson/Manu)
kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
                  "#D73202", "#272318", "#D3CDBF")

# --- load & label --------------------------------------------------------
df <- read_tsv(input_path, show_col_types = FALSE)
clans <- read_tsv(input_clan_path, show_col_types = FALSE)

required_cols <- c("stat_a_into_b", "stat_b_into_a")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "missing column(s) %s in %s -- available columns: %s",
    paste(missing_cols, collapse = ", "), input_path, paste(names(df), collapse = ", ")
  ))
}
df$combined_stat <- df$stat_a_into_b + df$stat_b_into_a

# Build a lookup of accession -> clan, then annotate each pair. Direct
# match on model_a/model_b against clans$accession, mirroring the
# minmax_boxplot.R pattern this script is based on -- no name/accession
# translation layer, since that script's working version does the same
# direct match.
clan_map <- setNames(clans$clan, clans$accession)
df$clan_a <- clan_map[df$model_a]
df$clan_b <- clan_map[df$model_b]

df$comparison <- ifelse(
  !is.na(df$clan_a) & !is.na(df$clan_b) & df$clan_a != "" & df$clan_a == df$clan_b,
  "within_clan",
  "everything_else"
)
df$comparison <- factor(df$comparison, levels = c("within_clan", "everything_else"))

if (all(is.na(df$clan_a)) && all(is.na(df$clan_b))) {
  warning("clan_a/clan_b are all NA -- model_a/model_b likely don't match clans$accession directly. ",
          "Check whether bitscoresumdifference's output uses family NAMEs rather than ACCESSIONs ",
          "(this differs from minmaxcmcompare's output if the two tools' population extraction diverged).")
}

# quick Wilcoxon (Mann-Whitney) test
wt <- wilcox.test(as.formula(paste(stat_column, "~ comparison")), data = df)
n_labels <- table(df$comparison)
subtitle_txt <- sprintf(
  "within_clan n=%d, everything_else n=%d  |  Wilcoxon p = %.3g",
  n_labels["within_clan"], n_labels["everything_else"], wt$p.value
)

# --- plot ------------------------------------------------------------
p <- ggplot(df, aes(x = comparison, y = .data[[stat_column]], fill = comparison)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.6, width = 0.5) +
  scale_fill_manual(values = c(within_clan = kea_palette[1],
                                everything_else = kea_palette[2])) +
  labs(
    x = NULL,
    y = "stat_a_into_b + stat_b_into_a",
    title = "Cross-recognition statistic: within-clan vs. everything else",
    subtitle = subtitle_txt,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p)
ggsave(output_path, p, width = 6, height = 5, dpi = 300)


# ---- modeling ----

df$within_clan <- as.integer(
  !is.na(df$clan_a) &
  !is.na(df$clan_b) &
  df$clan_a != "" &
  df$clan_a == df$clan_b
)

library(lme4)

m1 <- lmer(
  median_score_a_on_b ~ within_clan +
    (1 | model_a) +
    (1 | model_b),
  data = df
)

##summary(m1)
fixef(m1)

df$cross_sum <-
  df$median_score_a_on_b +
  df$median_score_b_on_a

m2 <- lmer(
  cross_sum ~ within_clan +
    (1 | model_a) +
    (1 | model_b),
  data = df
)

##summary(m2)
fixef(m2)

m3 <- lmer(
  stat_a_into_b ~ within_clan +
    (1 | model_a) +
    (1 | model_b),
  data = df
)

#summary(m3)
fixef(m3)

library(ggplot2)

ggplot(df,
       aes(median_score_b_on_b,
           median_score_a_on_b,
           color = factor(within_clan))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm")


df$cross_sum <-
    df$median_score_a_on_b +
    df$median_score_b_on_a

df$self_sum <-
    df$median_score_a_on_a +
    df$median_score_b_on_b

m4 <- lm(
  cross_sum ~ self_sum + within_clan,
  data = df
)
summary(m4)

m_cross <- lm(
    cross_sum ~ within_clan,
    data=df
)

summary(m_cross)

anova(m_cross, m4)



# library(readr)
# library(ggplot2)

# # --- config ------------------------------------------------------------
# input_path  <- "results/minmax_benchmark_5664258/raw_results.tsv"             # bitscoresumdifference/clan-benchmark output
# input_clan_path  <- "results/clan_population_5664258/population_clans.tsv"
# output_path <- "minmax_boxplot.png"

# # Manu "Kea" palette (github.com/G-Thomson/Manu)
# kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
#                   "#D73202", "#272318", "#D3CDBF")

# # --- load & label --------------------------------------------------------
# df <- read_tsv(input_path, show_col_types = FALSE)
# clans <- read_tsv(input_clan_path, show_col_types = FALSE)

# # Build a lookup of accession -> clan, then annotate each pair.
# clan_map <- setNames(clans$clan, clans$accession)
# df$clan_a <- clan_map[df$model_a]
# df$clan_b <- clan_map[df$model_b]

# df$comparison <- ifelse(
#   !is.na(df$clan_a) & !is.na(df$clan_b) & df$clan_a != "" & df$clan_a == df$clan_b,
#   "within_clan",
#   "everything_else"
# )

# df$comparison <- factor(df$comparison, levels = c("within_clan", "everything_else"))

# # quick Wilcoxon (Mann-Whitney) test
# wt <- wilcox.test(similarity_ratio ~ comparison, data = df)
# n_labels <- table(df$comparison)
# subtitle_txt <- sprintf(
#   "within_clan n=%d, everything_else n=%d  |  Wilcoxon p = %.3g",
#   n_labels["within_clan"], n_labels["everything_else"], wt$p.value
# )

# # --- plot ------------------------------------------------------------
# p <- ggplot(df, aes(x = comparison, y = similarity_ratio, fill = comparison)) +
#   geom_boxplot(outlier.shape = 21, outlier.alpha = 0.6, width = 0.5) +
#   scale_fill_manual(values = c(within_clan = kea_palette[1],
#                                 everything_else = kea_palette[2])) +
#   labs(
#     x = NULL,
#     y = "Similarity Ratio",
#     title = "Cross-recognition statistic: within-clan vs. everything else",
#     subtitle = subtitle_txt,
#     fill = NULL
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(legend.position = "none")

# print(p)
# ggsave(output_path, p, width = 6, height = 5, dpi = 300)