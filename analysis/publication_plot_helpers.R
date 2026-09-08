# Shared publication-style plotting helpers for benchmark analyses.
#
# Usage:
#   source("analysis/publication_plot_helpers.R")
#   p <- plot_clan_boxplot(df, stat_col = "stat",
#                          title = "Cross-recognition statistic",
#                          y_label = "KS D(a into b) + KS D(b into a)")
#
# Requirements: ggplot2

library(ggplot2)

publication_palette <- c(
  within_clan = "#6C803A",
  everything_else = "#7B5C34"
)

publication_labels <- c(
  within_clan = "Within clan",
  everything_else = "Other pairs"
)

theme_publication <- function(base_size = 12,
                             base_family = "",
                             legend_position = "none") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "#1f1f1f"),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.text = element_text(color = "#2d2d2d", size = rel(0.95)),
      axis.line = element_line(color = "#2d2d2d", linewidth = 0.6),
      axis.ticks = element_line(color = "#2d2d2d", linewidth = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e5e5e5", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(1.35), hjust = 0),
      plot.subtitle = element_text(color = "#4d4d4d", size = rel(0.9)),
      plot.caption = element_text(color = "#4d4d4d", size = rel(0.8)),
      legend.position = legend_position,
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text = element_text(size = rel(0.8)),
      legend.key.size = unit(0.55, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.spacing = unit(0.2, "cm"),
      legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(0.3, 0.4, 0.2, 0.2, "cm")
    )
}

prepare_clan_comparison <- function(df,
                                    stat_col,
                                    clans,
                                    accession_col = "accession",
                                    clan_col = "clan",
                                    model_a_col = "model_a",
                                    model_b_col = "model_b") {
  if (!stat_col %in% names(df)) {
    stop(sprintf("'%s' is not a column in the input data frame.", stat_col))
  }

  if (!is.data.frame(df)) {
    stop("df must be a data frame or tibble.")
  }

  df <- as.data.frame(df)

  if (missing(clans) || is.null(clans)) {
    stop("A clan table must be supplied to annotate within-clan vs. other-pair comparisons.")
  }

  clan_lookup <- stats::setNames(clans[[clan_col]], clans[[accession_col]])

  df$clan_a <- clan_lookup[df[[model_a_col]]]
  df$clan_b <- clan_lookup[df[[model_b_col]]]

  df$comparison <- ifelse(
    !is.na(df$clan_a) & !is.na(df$clan_b) &
      df$clan_a != "" & df$clan_a == df$clan_b,
    "within_clan",
    "everything_else"
  )

  df$comparison <- factor(df$comparison, levels = c("within_clan", "everything_else"))
  df
}

make_comparison_summary <- function(df, stat_col, comparison_col = "comparison") {
  if (!comparison_col %in% names(df)) {
    stop(sprintf("'%s' is not present in the data frame.", comparison_col))
  }

  if (!stat_col %in% names(df)) {
    stop(sprintf("'%s' is not present in the data frame.", stat_col))
  }

  wt <- wilcox.test(stats::as.formula(paste(stat_col, "~", comparison_col)), data = df)
  n_labels <- table(df[[comparison_col]])

  list(
    within_clan_n = as.integer(n_labels["within_clan"]),
    everything_else_n = as.integer(n_labels["everything_else"]),
    p_value = wt$p.value,
    p_text = sprintf("Wilcoxon p = %.3g", wt$p.value),
    n_text = sprintf("n = %d and %d", n_labels["within_clan"], n_labels["everything_else"]),
    summary_text = sprintf(
      "Within-clan n=%d; other pairs n=%d; Wilcoxon p = %.3g",
      n_labels["within_clan"],
      n_labels["everything_else"],
      wt$p.value
    )
  )
}

write_comparison_summary <- function(summary, output_path) {
  if (is.null(output_path) || output_path == "") {
    return(summary)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(
    sprintf("within_clan_n=%d", summary$within_clan_n),
    sprintf("everything_else_n=%d", summary$everything_else_n),
    sprintf("p_value=%.6g", summary$p_value),
    sprintf("p_text=%s", summary$p_text),
    sprintf("n_text=%s", summary$n_text),
    sprintf("summary_text=%s", summary$summary_text)
  )
  writeLines(lines, con = output_path)
  summary
}

plot_clan_boxplot <- function(df,
                              stat_col,
                              y_label,
                              title = NULL,
                              save_path = NULL,
                              stats_output_path = NULL,
                              width = 5.5,
                              height = 4.5,
                              dpi = 300,
                              box_width = 0.45,
                              y_limits = NULL,
                              fill_palette = publication_palette,
                              comparison_labels = publication_labels,
                              theme = theme_publication()) {
  if (!stat_col %in% names(df)) {
    stop(sprintf("'%s' was not found in the input data frame.", stat_col))
  }

  if (!"comparison" %in% names(df)) {
    stop("The input data frame must contain a 'comparison' column.")
  }

  stats_summary <- make_comparison_summary(df, stat_col = stat_col)
  if (!is.null(stats_output_path)) {
    write_comparison_summary(stats_summary, stats_output_path)
  }

  p <- ggplot(df, aes(x = comparison, y = .data[[stat_col]], fill = comparison)) +
    geom_boxplot(width = box_width, outlier.shape = 21, outlier.size = 1.5, outlier.alpha = 0.75) +
    scale_fill_manual(values = fill_palette, labels = comparison_labels, guide = "none") +
    scale_x_discrete(labels = comparison_labels) +
    labs(
      x = NULL,
      y = y_label,
      title = title,
      fill = NULL
    ) +
    theme

  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits)
  }

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }

  list(plot = p, summary = stats_summary)
}

plot_heatmap <- function(df, group_col, stat_col, reverse_stat_col = NULL, order_group_col = NULL, model_a_col = "model_a", model_b_col = "model_b", title = NULL, save_path = NULL, width = 6, height = 5, dpi = 300) {
  if (!is.data.frame(df)) {
    stop("df must be a data frame or tibble.")
  }

  required_cols <- c(model_a_col, model_b_col, group_col, stat_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")))
  }
  if (!is.null(reverse_stat_col) && !reverse_stat_col %in% names(df)) {
    stop(sprintf("'%s' was not found in the input data frame.", reverse_stat_col))
  }

  model_order <- unique(c(df[[model_a_col]], df[[model_b_col]]))
  boundary_positions <- numeric(0)
  if (!is.null(order_group_col)) {
    order_cols <- paste0(order_group_col, c("_a", "_b"))
    if (all(order_cols %in% names(df))) {
      model_groups <- rbind(
        data.frame(model = df[[model_a_col]], group = df[[order_cols[1]]]),
        data.frame(model = df[[model_b_col]], group = df[[order_cols[2]]])
      )
      model_groups <- model_groups[!duplicated(model_groups$model), ]
      model_groups$group[is.na(model_groups$group) | model_groups$group == ""] <- "Ungrouped"
      model_groups <- model_groups[order(model_groups$group, model_groups$model), ]
      model_order <- model_groups$model
      group_ends <- cumsum(table(factor(model_groups$group, levels = unique(model_groups$group))))
      boundary_positions <- as.numeric(group_ends)[-length(group_ends)] + 0.5
    }
  }
  y_boundary_positions <- length(model_order) + 1 - boundary_positions

  plot_df <- data.frame(
    model_x = df[[model_a_col]],
    model_y = df[[model_b_col]],
    value = df[[stat_col]],
    group = df[[group_col]]
  )
  if (!is.null(reverse_stat_col)) {
    plot_df <- rbind(plot_df, data.frame(
      model_x = df[[model_b_col]],
      model_y = df[[model_a_col]],
      value = df[[reverse_stat_col]],
      group = df[[group_col]]
    ))
  }

  p <- ggplot(plot_df, aes(
    x = factor(model_x, levels = model_order),
    y = factor(model_y, levels = rev(model_order)),
    fill = value,
    color = group
  )) +
    geom_tile(linewidth = 0.35) +
    scale_fill_viridis_c(option = "plasma", name = "Statistic") +
    scale_color_viridis_d(option = "turbo", name = group_col) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(hjust = 1, size = 7)
    ) +
    geom_vline(xintercept = boundary_positions, color = "white", linewidth = 0.8) +
    geom_hline(yintercept = y_boundary_positions, color = "white", linewidth = 0.8)

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(save_path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  }

  p
}