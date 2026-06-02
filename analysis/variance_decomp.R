#!/usr/bin/env Rscript
# Methods-note result: how is the spread of estimates across the multiverse
# apportioned among the analytic-choice dimensions? An ANOVA decomposition of the
# identified specification estimates by dimension, per estimand. Reads the engine's
# multiverse export (verification/js_multiverse.csv).
#
# Reports four complementary measures so the result does not hinge on a single,
# order-sensitive convention:
#   (i)   Type-I (sequential) SS share with estimator entered first  — "% sequential";
#   (ii)  Type-III (marginal) SS share, robust to ordering           — "% Type-III";
#   (iii) LMG / Shapley decomposition (averaged over all orderings)  — "% Shapley";
#   (iv)  the model R^2 from estimator alone                          — "estimator-only R^2";
# plus a residual = 100 - sum(shares) so each row sums to 100%, and the
# controls-only-vs-full spread contrast that the paper now quotes.

suppressMessages({library(readr); library(car); library(relaimpo)})
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
d <- read_csv("../verification/js_multiverse.csv", show_col_types = FALSE)
d <- d[d$identified == 1 & !is.na(d$b), ]

decomp <- function(sub) {
  sub$estimator <- factor(sub$estimator); sub$sample <- factor(sub$sample); sub$yearFE <- factor(sub$yearFE)
  full <- lm(b ~ estimator + experience + industry + region + health + yearFE + sample, data = sub)
  # Type-I sequential SS (estimator first)
  aI <- anova(full); ssI <- aI[["Sum Sq"]]; names(ssI) <- rownames(aI)
  totSS <- sum(ssI); shareI <- round(100 * ssI / totSS, 1)
  # Type-III marginal SS (order-insensitive)
  aIII <- car::Anova(full, type = 3); ssIII <- aIII[["Sum Sq"]]; names(ssIII) <- rownames(aIII)
  totIII <- sum(ssIII); shareIII <- round(100 * ssIII / totIII, 1)
  # LMG / Shapley decomposition (average marginal R^2 over all orderings)
  ri <- tryCatch(calc.relimp(full, type = "lmg", rela = FALSE), error = function(e) NULL)
  shareLMG <- if (!is.null(ri)) round(100 * ri@lmg, 1) else rep(NA_real_, 7)
  if (!is.null(ri)) names(shareLMG) <- names(ri@lmg)
  list(typeI = shareI, typeIII = shareIII, lmg = shareLMG,
       est_only_R2 = round(summary(lm(b ~ estimator, data = sub))$r.squared, 3),
       full_R2 = round(summary(full)$r.squared, 3),
       n = nrow(sub), range = range(sub$b), sd = round(sd(sub$b), 4))
}

# controls-only-vs-full multiverse contrast (the paper's new numerical sentence):
# fixing the estimator at pooled OLS, vary only controls + sample + yearFE, and
# compare the spread to the full multiverse including the estimator axis.
contrast <- function(sub) {
  full_lo <- min(sub$b); full_hi <- max(sub$b); full_iqr <- IQR(sub$b)
  pooled <- sub[sub$estimator == "pooled", ]
  if (nrow(pooled)) {
    co_lo <- min(pooled$b); co_hi <- max(pooled$b); co_iqr <- IQR(pooled$b)
  } else { co_lo <- NA; co_hi <- NA; co_iqr <- NA }
  list(full = c(full_lo, full_hi, full_iqr), controls_only = c(co_lo, co_hi, co_iqr))
}

for (f in c("union", "married", "educ")) {
  sub <- d[d$focal == f, ]
  r <- decomp(sub); cn <- contrast(sub)
  cat(sprintf("\n=== %s premium: %d identified specs, SD(estimate)=%.3f, range=[%.3f, %.3f] ===\n", f, r$n, r$sd, r$range[1], r$range[2]))
  cat("Type-I sequential SS shares (estimator first; %):\n"); print(r$typeI[names(r$typeI) != "Residuals"])
  cat(sprintf("  Residual (within-cells noise): %.1f%%\n", r$typeI[["Residuals"]]))
  cat("Type-III marginal SS shares (order-insensitive; %):\n"); print(r$typeIII[!(names(r$typeIII) %in% c("Residuals", "(Intercept)"))])
  cat("LMG / Shapley shares (averaged over all orderings; %):\n"); print(r$lmg)
  cat(sprintf("  estimator-only R^2 = %.3f; full-model R^2 = %.3f\n", r$est_only_R2, r$full_R2))
  cat(sprintf("Controls-only spread (estimator fixed at pooled OLS): [%.3f, %.3f], IQR %.3f\n", cn$controls_only[1], cn$controls_only[2], cn$controls_only[3]))
  cat(sprintf("Full-multiverse  spread (all estimators):              [%.3f, %.3f], IQR %.3f\n", cn$full[1], cn$full[2], cn$full[3]))
}
