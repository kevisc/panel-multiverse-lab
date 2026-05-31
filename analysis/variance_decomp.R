#!/usr/bin/env Rscript
# Methods-note result: how is the spread of estimates across the multiverse
# apportioned among the analytic-choice dimensions? An ANOVA decomposition of the
# identified specification estimates by dimension, per estimand. Reads the engine's
# multiverse export (verification/js_multiverse.csv).
suppressMessages({library(readr)})
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
d <- read_csv("../verification/js_multiverse.csv", show_col_types = FALSE)
d <- d[d$identified == 1 & !is.na(d$b), ]

decomp <- function(sub) {
  sub$estimator <- factor(sub$estimator); sub$sample <- factor(sub$sample)
  sub$yearFE <- factor(sub$yearFE)
  # main-effects ANOVA of the point estimate on the choice dimensions
  m <- lm(b ~ estimator + experience + industry + region + health + yearFE + sample, data = sub)
  a <- anova(m); ss <- a[["Sum Sq"]]; names(ss) <- rownames(a)
  tot <- sum(ss)
  share <- round(100 * ss / tot, 1)
  list(share = share, R2 = summary(m)$r.squared, n = nrow(sub),
       range = range(sub$b), sd = sd(sub$b),
       est_only_R2 = summary(lm(b ~ estimator, data = sub))$r.squared)
}

for (f in c("union", "married", "educ")) {
  sub <- d[d$focal == f, ]
  r <- decomp(sub)
  cat(sprintf("\n=== %s premium: %d identified specs, SD(estimate)=%.3f, range=[%.3f, %.3f] ===\n",
              f, r$n, r$sd, r$range[1], r$range[2]))
  cat("Share of estimate variance (ANOVA Sum-of-Squares, %):\n")
  print(r$share[names(r$share) != "Residuals"])
  cat(sprintf("  estimator alone explains R^2 = %.3f of the variance in estimates\n", r$est_only_R2))
}
