#!/usr/bin/env Rscript
# Builds the specification-curve ("universe") figure for the paper from the
# multiverse exported by the JS engine (figures/multiverse_union.csv, produced
# by enumerating the union-premium multiverse in engine.js). Publication theme.
suppressMessages({library(ggplot2); library(patchwork); library(readr)})
setwd(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))

d <- read_csv("multiverse_union.csv", show_col_types = FALSE)
d$rank <- seq_len(nrow(d))
d$Significant <- factor(ifelse(d$sig == 1, "95% CI excludes 0", "not significant"),
                        levels = c("95% CI excludes 0", "not significant"))
med <- median(d$b)
cols <- c("95% CI excludes 0" = "#1f6f3d", "not significant" = "#9aa0a6")

# ---- top panel: sorted estimates with confidence intervals ----
top <- ggplot(d, aes(rank, b, color = Significant)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = med, linetype = "dashed", color = "#b22222", linewidth = 0.5) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.3, alpha = 0.6) +
  geom_point(size = 0.7) +
  annotate("text", x = nrow(d), y = med, label = sprintf("median = %.3f", med),
           hjust = 1, vjust = -0.6, size = 3, color = "#b22222") +
  scale_color_manual(values = cols, name = NULL) +
  labs(y = "Union wage premium\n(log points)") +
  theme_minimal(base_size = 11) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        legend.position = "top", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

# ---- bottom panel: analytic-choice dashboard ----
est_lab <- c(pooled = "Pooled OLS", between = "Between", random = "Random effects",
             within = "Fixed effects", fd = "First differences")
rowsdef <- rbind(
  data.frame(label = est_lab[c("pooled","between","random","within","fd")],
             on = sapply(c("pooled","between","random","within","fd"),
                         function(e) list(d$estimator == e))),
  NULL)
# Build a long data frame of active cells
mk <- function(label, active) data.frame(rank = d$rank, label = label, active = active, sig = d$sig)
rows_list <- list(
  mk("Pooled OLS",        d$estimator == "pooled"),
  mk("Between",           d$estimator == "between"),
  mk("Random effects",    d$estimator == "random"),
  mk("Fixed effects",     d$estimator == "within"),
  mk("First differences", d$estimator == "fd"),
  mk("experience (+sq.)", d$experience == 1),
  mk("industry",          d$industry == 1),
  mk("region",            d$region == 1),
  mk("health",            d$health == 1),
  mk("year effects",      d$yearFE == 1)
)
dash <- do.call(rbind, rows_list)
lev <- rev(c("Pooled OLS","Between","Random effects","Fixed effects","First differences",
             "experience (+sq.)","industry","region","health","year effects"))
dash$label <- factor(dash$label, levels = lev)
dash <- dash[dash$active, ]
dash$Significant <- factor(ifelse(dash$sig == 1, "95% CI excludes 0", "not significant"),
                           levels = c("95% CI excludes 0", "not significant"))

bottom <- ggplot(dash, aes(rank, label, color = Significant)) +
  geom_point(size = 0.7, show.legend = FALSE) +
  scale_color_manual(values = cols) +
  labs(x = "specifications, ordered by estimate", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 9))

fig <- top / bottom + plot_layout(heights = c(1.5, 1.6))

ggsave("spec-curve-union.png", fig, width = 7.2, height = 6.0, dpi = 300, bg = "white")
cat("wrote figures/spec-curve-union.png\n")
