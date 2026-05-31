#!/usr/bin/env Rscript
# Figure 2: the estimator "spectrum". For each estimand, the distribution of the focal
# coefficient across all identified specifications, by estimator. Shows the monotone
# within < between ordering for the union/marriage premiums and the absent fixed-effects /
# first-difference estimates for the (time-invariant) schooling return.
suppressMessages({library(ggplot2); library(readr)})
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

d <- read_csv("../verification/js_multiverse.csv", show_col_types = FALSE)
d <- d[d$identified == 1 & !is.na(d$b), ]
d$sig <- (d$b - 1.96*d$se > 0) | (d$b + 1.96*d$se < 0)
d$Significant <- factor(ifelse(d$sig, "95% CI excludes 0", "not significant"),
                        levels = c("95% CI excludes 0", "not significant"))
est_lab <- c(pooled="Pooled\nOLS", between="Between", random="Random\neffects",
             within="Fixed\neffects", fd="First\ndiff.")
d$estimator <- factor(est_lab[d$estimator], levels = est_lab[c("pooled","between","random","within","fd")])
foc_lab <- c(union="Union premium", married="Marriage premium", educ="Return to schooling")
d$focal <- factor(foc_lab[d$focal], levels = foc_lab[c("union","married","educ")])
cols <- c("95% CI excludes 0" = "#1f6f3d", "not significant" = "#9aa0a6")

p <- ggplot(d, aes(estimator, b, color = Significant)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_jitter(width = 0.18, height = 0, size = 0.6, alpha = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, color = "black",
               linewidth = 0.3, fatten = 0) +
  facet_wrap(~ focal, nrow = 1) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = NULL, y = "Focal coefficient (log points)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(face = "bold"))

ggsave("estimator-spectrum.png", p, width = 7.4, height = 3.6, dpi = 300, bg = "white")
cat("wrote figures/estimator-spectrum.png\n")
# quick console summary of the missing schooling FE/FD
cat("Return-to-schooling specs by estimator:\n")
print(table(d$focal, d$estimator)["Return to schooling", ])
