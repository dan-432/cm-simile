# cluster_candidates.R
# Usage: Rscript cluster_candidates.R cmsimile_pairwise.tsv candidate_clustering
# Writes: <prefix>_dendrogram.pdf and <prefix>_heatmap.pdf

args     <- commandArgs(trailingOnly = TRUE)
in_path  <- if (length(args) >= 1) args[1] else "cmsimile_pairwise.tsv"
prefix   <- if (length(args) >= 2) args[2] else "candidate_clustering"
helper   <- "analysis/publication_plot_helpers.R"
linkage  <- "average"   # alternatives: "complete", "ward.D2"

message("Working directory: ", getwd())
if (!file.exists(in_path)) stop("Input file not found: ", in_path)
if (!file.exists(helper))  stop("Helper script not found: ", helper)

source(helper)   # also loads ggplot2 and defines theme_publication
df <- read.delim(in_path, stringsAsFactors = FALSE)
required_cols <- c("model_a", "model_b", "ks_d_a_into_b", "ks_d_b_into_a")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

short_name <- function(x) {
  x <- sub("^model_[0-9]+_cluster_([0-9]+)_.*$", "\\1", x)
  x <- sub("^(RF[0-9]+)_.*$", "\\1", x)
  sub("_1\\.cacofold$", "", x)
}
df$a <- short_name(df$model_a)
df$b <- short_name(df$model_b)
df$d_sum <- df$ks_d_a_into_b + df$ks_d_b_into_a   # similarity, 0 to 2

models <- sort(unique(c(df$a, df$b)))
sim <- matrix(2, nrow = length(models), ncol = length(models),
              dimnames = list(models, models))
sim[cbind(df$a, df$b)] <- df$d_sum
sim[cbind(df$b, df$a)] <- df$d_sum

# Keep the clustering symmetric, but show the directional KS D value in
# each heatmap cell. The mirrored cell therefore shows the reverse direction.
directional <- matrix(NA_real_, nrow = length(models), ncol = length(models),
                      dimnames = list(models, models))
directional[cbind(df$a, df$b)] <- df$ks_d_a_into_b
directional[cbind(df$b, df$a)] <- df$ks_d_b_into_a

hc <- hclust(as.dist(2-sim), method = linkage)

dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)

# ---- Dendrogram (no title) ----
dendro_path <- paste0(prefix, "_dendrogram.pdf")
pdf(dendro_path, width = 5.5, height = 4)
tryCatch({
  par(mar = c(1, 4.5, 1, 1))
  plot(hc, hang = -1, xlab = "", sub = "", main = "",
       ylab = "Distance: 2 - (KS D a->b + KS D b->a)")
}, finally = dev.off())

# ---- Heatmap in dendrogram order (no title) ----
ord <- hc$labels[hc$order]
plot_df <- expand.grid(x = ord, y = ord, stringsAsFactors = FALSE)
plot_df$value <- directional[cbind(plot_df$x, plot_df$y)]
plot_df$value[plot_df$x == plot_df$y] <- NA

p <- ggplot(plot_df, aes(x = factor(x, levels = ord),
                         y = factor(y, levels = rev(ord)),
                         fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1),
                       na.value = "grey90", name = "Directional KS D") +
  labs(x = NULL, y = NULL) +
  theme_publication(legend_position = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

heat_path <- paste0(prefix, "_heatmap.pdf")
ggsave(heat_path, plot = p, width = 6, height = 4, device = "pdf", bg = "white")

message("Wrote ", normalizePath(dendro_path))
message("Wrote ", normalizePath(heat_path))