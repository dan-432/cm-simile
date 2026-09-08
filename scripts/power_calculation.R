#!/usr/bin/env Rscript
# Power calculation for cm-simile's comparisons, based on sample size and
# variance -- built on base R's power.t.test() (stats package, zero extra
# dependencies), which is both directly applicable and the standard
# practical approximation for Mann-Whitney/Wilcoxon power via the test's
# asymptotic relative efficiency (ARE) relative to the t-test.
#
# Why a t-test-based calculation is a reasonable stand-in for Mann-
# Whitney power, not a mismatch: under (approximately) normal data, the
# Wilcoxon rank-sum test's ARE relative to the t-test is 3/pi ~= 0.955 --
# it needs about 1/0.955 ~= 4.7% more samples to match the t-test's power
# for the same effect size. So a t-test power calculation, with sample
# size scaled by that factor, gives a solid approximation for Mann-
# Whitney power without needing an exact (and much less "easy") rank-
# based power formula. This is the standard applied-stats shortcut, not
# a hack -- it degrades gracefully as data becomes non-normal (Mann-
# Whitney's actual ARE advantage over the t-test only grows for heavier-
# tailed distributions, so this approximation is, if anything,
# conservative).
#
# Two modes, both via the same underlying calculation:
#   1. Given n (per group) + effect size (mean difference, SD) -> power
#   2. Given target power + effect size -> required n (per group)
# Leave exactly one of `n` or `power` as NA below to choose which way to
# solve.
#
# Requires: base R only (stats package, always available).

# --- config --------------------------------------------------------------
# Effect size inputs -- estimate these from your actual data, e.g. the
# real-pool-vs-null mean/SD difference for a specific pair or averaged
# across several, from real_matrix.tsv / reservoir_matrix.tsv.
mean_diff   <- 10       # expected difference in means between real and null scores
pooled_sd   <- 15       # pooled (or either sample's, if roughly equal) standard deviation

# Sample size / power -- set exactly one to NA.
n     <- 100            # per-group sample size (e.g. -npool or -resample-size); NA to solve for n
power <- NA             # target power (e.g. 0.8); NA to solve for power given n above

alpha <- 0.05            # two-sided significance threshold
mannwhitney_are <- 3 / pi # asymptotic relative efficiency of Wilcoxon rank-sum vs. t-test under normality

# --- calculation -----------------------------------------------------

if (is.na(n) == is.na(power)) {
  stop("set exactly one of `n` or `power` to NA -- the other is the target you're solving for")
}

d <- mean_diff / pooled_sd  # Cohen's d

if (is.na(n)) {
  # Solve for required n (t-test), then inflate by 1/ARE for the
  # Mann-Whitney-equivalent requirement -- Wilcoxon needs somewhat more
  # samples than the t-test to hit the same power for the same effect.
  t_result <- power.t.test(delta = mean_diff, sd = pooled_sd, sig.level = alpha, power = power, type = "two.sample")
  n_t  <- t_result$n
  n_mw <- n_t / mannwhitney_are
  cat(sprintf("Effect size (Cohen's d): %.3f\n", d))
  cat(sprintf("Required n per group for power=%.2f, alpha=%.3f:\n", power, alpha))
  cat(sprintf("  t-test:        n = %.1f\n", n_t))
  cat(sprintf("  Mann-Whitney:  n = %.1f  (t-test n / ARE, ARE=%.4f)\n", n_mw, mannwhitney_are))
} else {
  # Solve for power at the given n (t-test), then compute the
  # Mann-Whitney-equivalent power by scaling n down by ARE first (a
  # smaller effective sample size for the same nominal n, since
  # Wilcoxon is slightly less efficient than the t-test under normality).
  t_result <- power.t.test(n = n, delta = mean_diff, sd = pooled_sd, sig.level = alpha, type = "two.sample")
  mw_result <- power.t.test(n = n * mannwhitney_are, delta = mean_diff, sd = pooled_sd, sig.level = alpha, type = "two.sample")
  cat(sprintf("Effect size (Cohen's d): %.3f\n", d))
  cat(sprintf("Power at n=%d per group, alpha=%.3f:\n", n, alpha))
  cat(sprintf("  t-test:        power = %.4f\n", t_result$power))
  cat(sprintf("  Mann-Whitney:  power = %.4f  (ARE-adjusted, ARE=%.4f)\n", mw_result$power, mannwhitney_are))
}

cat("\nNote: this is the t-test/Wilcoxon rank-sum power calculation for a\n")
cat("SINGLE pair comparison at the given alpha. If you're running this\n")
cat("across many pairs and using an FDR/Bonferroni-adjusted threshold\n")
cat("rather than raw alpha, substitute your adjusted alpha above --\n")
cat("power drops as the effective significance threshold tightens.\n")

cat("\nThis does NOT give power for the Kolmogorov-Smirnov test (ks_d_*\n")
cat("columns from ksmwcompare) -- KS tests general distributional\n")
cat("difference, not a location shift, and doesn't have as simple a\n")
cat("closed-form power calculation from just mean difference + SD.\n")
cat("If you need KS power specifically, simulation (draw many synthetic\n")
cat("real/null pairs at your candidate n and effect size, run TwoSampleKS\n")
cat("on each, count how often p < alpha) is the standard approach --\n")
cat("ask if you'd like a script for that.\n")
