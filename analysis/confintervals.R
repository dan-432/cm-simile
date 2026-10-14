#!/usr/bin/env Rscript
# Bootstrap-over-families confidence intervals for the clan benchmark AUC

library(readr)
source(file.path("analysis", "publication_plot_helpers.R"))
source(file.path("analysis", "aoc_stats_per_approach.R"))   # defines `approaches`, as in compute_auc.R

n_boot <- 5000
set.seed(1)

clans <- read_tsv("results/clan_population_5687123/population_clans.tsv",
                  show_col_types = FALSE)

# Oriented statistic and clan label for every pair, per approach
prep <- lapply(approaches, function(a) {
  df <- read_tsv(a$path, show_col_types = FALSE)
  df$stat <- a$stat(df)
  prepare_clan_comparison(df, stat_col = "stat", clans = clans)
})
names(prep) <- vapply(approaches, function(a) a$name, character(1))

families <- sort(unique(c(prep[[1]]$model_a, prep[[1]]$model_b)))
n_fam <- length(families)
for (p in prep) stopifnot(setequal(c(p$model_a, p$model_b), families))

# Family resampling shared by all approaches: counts[b, f] = times family f is drawn
counts <- t(rmultinom(n_boot, size = n_fam, prob = rep(1 / n_fam, n_fam)))

boot_auc <- function(d) {
  within <- d$comparison == "within_clan"
  s_w <- d$stat[within]
  s_o <- d$stat[!within]
  # M[i, j] = 1 if within-clan pair i outscores other pair j, 0.5 for ties
  M <- outer(s_w, s_o, function(x, y) (x > y) + 0.5 * (x == y))
  ia <- match(d$model_a, families)
  ib <- match(d$model_b, families)

  reps <- vapply(seq_len(n_boot), function(b) {
    w  <- counts[b, ia] * counts[b, ib]   # a pair is kept once per copy of each family
    ww <- w[within]
    wo <- w[!within]
    if (sum(ww) == 0 || sum(wo) == 0) return(NA_real_)
    drop(ww %*% M %*% wo) / (sum(ww) * sum(wo))
  }, numeric(1))

  list(est = mean(M), reps = reps)
}

res <- lapply(prep, boot_auc)
ci  <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)

for (nm in names(res)) {
  r <- res[[nm]]
  q <- ci(r$reps)
  cat(sprintf("%s\n  AUC = %.3f (95%% CI %.3f to %.3f); %d of %d resamples usable\n\n",
              nm, r$est, q[1], q[2], sum(!is.na(r$reps)), n_boot))
}

nms <- names(res)
for (i in seq_len(length(nms) - 1)) {
  for (j in (i + 1):length(nms)) {
    q <- ci(res[[i]]$reps - res[[j]]$reps)
    cat(sprintf("%s minus %s\n  difference = %.3f (95%% CI %.3f to %.3f)\n\n",
                nms[i], nms[j], res[[i]]$est - res[[j]]$est, q[1], q[2]))
  }
}