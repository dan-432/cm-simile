#!/usr/bin/env Rscript
# Boxplot of stat_a_into_b, within-clan pairs vs. everything else.
#
# Requires: readr, ggplot2
#   install.packages(c("readr", "ggplot2"))
#   -- or, if you keep this alongside the cmsimile conda env --
#   conda install -c conda-forge r-base r-readr r-ggplot2

library(readr)
source(file.path("analysis", "publication_plot_helpers.R"))

# --- config ------------------------------------------------------------
input_path  <- "results/ksmwcompare_benchmark_5740367/raw_results.tsv"
input_clan_path  <- "results/clan_population_5740367/population_clans.tsv"
output_path <- file.path("analysis", "ksmwcompare_boxplot.svg")
stats_output_path <- file.path("analysis", "ksmwcompare_boxplot_stats.txt")

# --- load & label --------------------------------------------------------
df <- read_tsv(input_path, show_col_types = FALSE)
clans <- read_tsv(input_clan_path, show_col_types = FALSE)

df$stat <- df$ks_d_a_into_b + df$ks_d_b_into_a

df <- prepare_clan_comparison(
  df,
  stat_col = "stat",
  clans = clans,
  accession_col = "accession",
  clan_col = "clan",
  model_a_col = "model_a",
  model_b_col = "model_b"
)

# --- plot ------------------------------------------------------------
result <- plot_clan_boxplot(
  df,
  stat_col = "stat",
  y_label = "KS D(a into b) + KS D(b into a)",
  title = NULL,
  save_path = output_path,
  stats_output_path = stats_output_path,
  width = 5.5,
  height = 4.5,
  box_width = 0.45,
  y_limits = c(0, 2),
  dpi = 300
)

print(result$plot)
print(result$summary$summary_text)