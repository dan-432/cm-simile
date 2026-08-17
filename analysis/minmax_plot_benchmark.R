#!/usr/bin/env Rscript
# Boxplot of stat_a_into_b, within-clan pairs vs. everything else.
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
library(ggplot2)

# --- config ------------------------------------------------------------
input_path  <- "results/minmax_benchmark_5664258/raw_results.tsv"             # bitscoresumdifference/clan-benchmark output
input_clan_path  <- "results/clan_population_5664258/population_clans.tsv"
output_path <- "minmax_boxplot.png"

# Manu "Kea" palette (github.com/G-Thomson/Manu)
kea_palette <- c("#6C803A", "#7B5C34", "#AB7C47", "#CCAE42",
                  "#D73202", "#272318", "#D3CDBF")

# --- load & label --------------------------------------------------------
df <- read_tsv(input_path, show_col_types = FALSE)
clans <- read_tsv(input_clan_path, show_col_types = FALSE)

# Build a lookup of accession -> clan, then annotate each pair.
clan_map <- setNames(clans$clan, clans$accession)
df$clan_a <- clan_map[df$model_a]
df$clan_b <- clan_map[df$model_b]

df$comparison <- ifelse(
  !is.na(df$clan_a) & !is.na(df$clan_b) & df$clan_a != "" & df$clan_a == df$clan_b,
  "within_clan",
  "everything_else"
)

df$comparison <- factor(df$comparison, levels = c("within_clan", "everything_else"))

# quick Wilcoxon (Mann-Whitney) test
wt <- wilcox.test(similarity_ratio ~ comparison, data = df)
n_labels <- table(df$comparison)
subtitle_txt <- sprintf(
  "within_clan n=%d, everything_else n=%d  |  Wilcoxon p = %.3g",
  n_labels["within_clan"], n_labels["everything_else"], wt$p.value
)

# --- plot ------------------------------------------------------------
p <- ggplot(df, aes(x = comparison, y = similarity_ratio, fill = comparison)) +
  geom_boxplot(outlier.shape = 21, outlier.alpha = 0.6, width = 0.5) +
  scale_fill_manual(values = c(within_clan = kea_palette[1],
                                everything_else = kea_palette[2])) +
  labs(
    x = NULL,
    y = "Similarity Ratio",
    title = "Cross-recognition statistic: within-clan vs. everything else",
    subtitle = subtitle_txt,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p)
ggsave(output_path, p, width = 6, height = 5, dpi = 300)