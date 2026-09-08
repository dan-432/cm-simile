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
input_path  <- "results/ksmw-rf00050-allclust.tsv"
output_path <- file.path("analysis", "ksmwrf0050split_heatmap.svg")
stats_output_path <- file.path("analysis", "ksmwrf0050split_heatmap_stats.txt")

df <- read_tsv(input_path, show_col_types = FALSE)

df$comparison <- "none"

# --- plot ------------------------------------------------------------
result <- plot_heatmap(
  df,
  group_col = "comparison",
  stat_col = "ks_d_a_into_b",
  reverse_stat_col = "ks_d_b_into_a",
  # order_group_col = "clan",
  model_a_col = "model_a",
  model_b_col = "model_b",
  save_path = output_path,
  width = 9,
  height = 9,
  dpi = 300
)