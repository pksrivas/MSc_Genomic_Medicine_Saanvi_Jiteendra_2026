# Figure 5. Per gene change in editing ratio after ADAR knockdown, endothelial2.
#
# For every gene, the change in editing after knockdown is worked out as the
# gene's mean editing ratio in the knockdown samples minus its mean in the
# scrambled controls. A negative value means the gene was edited less after
# knockdown, which is what ADAR1 knockdown is expected to produce.
#
# The figure shows the whole distribution of those per gene values rather than
# summarising them in a table. A table of mean, median, interquartile range and
# percentage negative is five numbers standing in for a shape, and the shape is
# what the argument rests on. Each row here is one gene set: a density curve on
# top, a box underneath showing median and interquartile range, and a dashed
# line at zero so the proportion of genes below zero can be seen directly.
#
# FOUR GENE SETS, TWO COMPARISONS
#
# Rows are, from the top, all tested genes, then genes with at least one
# GLMM-significant site, then at least one Fisher-significant site, then genes
# with a single site called by both tests. The two panels are the two
# knockdowns.
#
# THE SELECTION EFFECT, WHICH IS THE MAIN THING TO UNDERSTAND HERE
#
# The three significance-selected rows look tidier than the top row. They sit
# further from zero and are almost entirely negative. That is largely built in
# rather than discovered: a site becomes significant partly because its effect
# is large, so selecting on significance and then observing large effects is
# close to circular.
#
# The row that carries the actual directional claim is the top one, all tested
# genes, because nothing was selected for. This caveat is printed in the
# figure's subtitle rather than left to the caption.
#
# WHY BOTH PANELS SHARE ONE HORIZONTAL SCALE
#
# The finding is that siADAR1 shifts left and siADAR2 barely moves. If each
# panel were free to set its own scale, siADAR2's small random variation would
# be stretched to fill its panel and would look just as pronounced as siADAR1's
# real shift. A shared scale keeps the comparison honest.
#
# THREE CHOICES ABOUT DRAWING SMALL GROUPS
#
# The gene sets range from 725 genes down to 4. Densities are drawn at a fixed
# width rather than scaled by group size, because scaling by size would make the
# small sets invisible. The group size is printed next to each row instead, so
# nothing is hidden by that choice.
#
# Sets with fewer than 20 genes are drawn as individual points rather than a
# density curve. A smooth curve through four values shows the shape of the
# smoothing calculation, not the shape of the data.
#
# Sets with no genes at all get an explicit label saying so, rather than an
# empty row that a reader might take for missing or failed data.
#
# WHY THE BOTTOM ROW SAYS "BOTH, SAME SITE"
#
# This row counts genes with at least one site called significant by both tests
# at once, which gives 168 genes. That is not the same as the 172 in the Euler
# diagram of Figure 1, which counts genes having at least one GLMM-significant
# site and at least one Fisher-significant site without requiring them to be the
# same site. Both counts are correct and they answer different questions, so
# this row is labelled by its definition rather than loosely as "agreement".
#
# COLOUR
#
# One colour per gene set, checked for readability and for colour vision
# deficiency. Teal is used for the bottom row rather than purple, even though
# purple separates slightly better, because purple already means the Wilcoxon
# test in Figure 6 and one colour should not stand for two things across a
# thesis. The top row is drawn in grey because it is context rather than a
# finding. Every row is also labelled in words, so nothing depends on colour
# alone.
#
# Input:  gene_table_siADAR1.txt, gene_table_siADAR2.txt
# Output: delta_distributions_endo2.png

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
O <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(O)

CTX <- "#8a897f"; GLMM_C <- "#2a78d6"; FISH_C <- "#eb6834"; BOTH_C <- "#0f9488"
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; MUTED <- "#52514e"; GRIDS <- "#e1e0d9"; TXT <- "#898781"

SETS <- c("All tested", "GLMM-sig", "Fisher-sig", "Both, same site")
COLS <- c(`All tested` = CTX, `GLMM-sig` = GLMM_C, `Fisher-sig` = FISH_C,
          `Both, same site` = BOTH_C)

# ---------------------------------------------------------------------------
# Assemble one long table: one row per gene, per gene set, per comparison
# ---------------------------------------------------------------------------

# A gene can belong to several sets at once, so it appears once per set it
# qualifies for. That is intended: the sets are nested views of the same genes,
# not a partition.
gather <- function(case) {
  g <- fread(sprintf("gene_table_%s.txt", case))
  # One membership test per set, in the same order as SETS above. The first is
  # every gene with a real symbol, which is the unselected reference set.
  sel <- list(g$SYMBOL != "", g$n_GLMM_sig >= 1, g$n_Fisher_sig >= 1, g$n_agree >= 1)
  rbindlist(lapply(seq_along(SETS), function(i)
    if (!any(sel[[i]])) data.table(comp = case, set = SETS[i], delta = numeric(0))
    else data.table(comp = case, set = SETS[i], delta = g$mean_delta[sel[[i]]])))
}
d <- rbindlist(lapply(c("siADAR1", "siADAR2"), gather))
d[, comp := factor(fifelse(comp == "siADAR1", "siADAR1 vs scr", "siADAR2 vs scr"),
                   levels = c("siADAR1 vs scr", "siADAR2 vs scr"))]
# reversed so "All tested" sits at the top of each panel
d[, set := factor(set, levels = rev(SETS))]

# Build every combination of comparison and gene set, including ones with no
# genes. Joining the real data onto this grid keeps empty sets as visible rows
# with a count of zero, instead of dropping them out of the figure entirely.
grid <- CJ(comp = levels(d$comp), set = levels(d$set), unique = TRUE)
grid[, comp := factor(comp, levels = levels(d$comp))][, set := factor(set, levels = levels(d$set))]
grid[, yi := as.integer(set)]

stat <- d[, .(n = .N, pneg = 100 * mean(delta < 0), mu = mean(delta),
              med = median(delta), q1 = quantile(delta, .25), q3 = quantile(delta, .75),
              lo = quantile(delta, .05), hi = quantile(delta, .95)), by = .(comp, set)]
stat <- merge(grid, stat, by = c("comp", "set"), all.x = TRUE)
stat[is.na(n), n := 0L]

# ---------------------------------------------------------------------------
# Density curves, drawn only for sets with enough genes to justify one
# ---------------------------------------------------------------------------

# Below this many genes, a density curve would describe the smoothing rather
# than the data, so those sets are drawn as individual points further down.
BIG <- 20L
viol <- rbindlist(lapply(split(d, list(d$comp, d$set), drop = TRUE), function(s) {
  if (nrow(s) < BIG) return(NULL)
  # density() normally extends past the smallest and largest observed values, so
  # the curve would imply genes exist beyond the range actually measured. Fixing
  # the limits to the observed minimum and maximum stops the curve claiming
  # data that is not there, and removes a flat hairline at the baseline that
  # otherwise reads as a stray horizontal rule.
  k <- density(s$delta, adjust = 0.9, from = min(s$delta), to = max(s$delta))
  yi <- as.integer(s$set[1])
  data.table(comp = s$comp[1], set = s$set[1],
             x = c(k$x, rev(k$x)),
             y = c(yi + 0.07 + 0.36 * k$y / max(k$y), rep(yi + 0.07, length(k$x))))
}))
pts <- merge(d, stat[, .(comp, set, n)], by = c("comp", "set"))[n < BIG]
pts[, yi := as.integer(set)]
set.seed(1); pts[, yj := yi - 0.10 + runif(.N, -0.05, 0.05)]

# ---------------------------------------------------------------------------
# A reserved strip on the right for printed numbers
# ---------------------------------------------------------------------------

# The group size and the percentage of genes below zero are printed beside each
# row rather than left to the reader to estimate. They sit in a strip on the
# right, separated from the plotting area by a vertical rule, and the axis
# gridlines are stopped short of it so nothing runs underneath the numbers.
X_LO <- -0.63; X_RULE <- 0.50; X_N <- 0.60; X_P <- 0.79; X_HI <- 0.88
stat[, `:=`(lab_n = format(n, big.mark = ","),
            lab_p = fifelse(n == 0, "--", sprintf("%.1f%%", pneg)))]
hdr <- data.table(comp = factor(levels(d$comp)[1], levels = levels(d$comp)),
                  x = c(X_N, X_P), y = 4.62, t = c("n", "% neg"))
empty <- stat[n == 0]

p <- ggplot() +
  annotate("rect", xmin = X_RULE, xmax = X_HI, ymin = -Inf, ymax = Inf,
           fill = SURF, colour = NA) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.55) +
  annotate("segment", x = X_RULE, xend = X_RULE, y = 0.45, yend = 4.45,
           colour = GRIDS, linewidth = 0.5) +
  # 5th-95th whisker, IQR box, median tick, mean diamond
  geom_segment(data = stat[n >= BIG], aes(x = lo, xend = hi, y = yi - 0.19, yend = yi - 0.19),
               colour = MUTED, linewidth = 0.4) +
  geom_rect(data = stat[n >= BIG],
            aes(xmin = q1, xmax = q3, ymin = yi - 0.30, ymax = yi - 0.08, fill = set),
            colour = SURF, linewidth = 0.6) +
  geom_polygon(data = viol, aes(x, y, group = interaction(comp, set), fill = set),
               alpha = 0.55, colour = NA) +
  # median is the white tick only. An additional median circle was drawn here
  # and removed: mean and median differ by <0.007 in every set, so the two
  # markers overlapped into a single unreadable blob.
  geom_segment(data = stat[n >= BIG], aes(x = med, xend = med, y = yi - 0.31, yend = yi - 0.07),
               colour = SURF, linewidth = 1.6) +
  geom_point(data = stat[n >= BIG], aes(mu, yi - 0.19), shape = 23, size = 2.7,
             fill = INK, colour = SURF, stroke = 0.6) +
  # small sets: raw points only
  geom_point(data = pts, aes(delta, yj, fill = set), shape = 21, size = 2.6,
             colour = SURF, stroke = 0.6) +
  # offset well clear of the zero rule, which it previously sat on top of
  geom_text(data = empty, aes(x = -0.34, y = yi - 0.10, label = "no genes in set"),
            hjust = 0.5, size = 3.05, colour = TXT, fontface = "italic") +
  # gutter
  geom_text(data = stat, aes(X_N, yi - 0.08, label = lab_n), hjust = 1,
            size = 3.3, colour = INK) +
  geom_text(data = stat, aes(X_P, yi - 0.08, label = lab_p), hjust = 1,
            size = 3.3, colour = INK) +
  geom_text(data = hdr, aes(x, y, label = t), hjust = 1, size = 3.05,
            colour = MUTED, fontface = "bold") +
  facet_wrap(~ comp, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_y_continuous(breaks = 1:4, labels = levels(d$set),
                     limits = c(0.42, 4.78), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-0.6, 0.4, 0.2),
                     labels = scales::number_format(accuracy = 0.1),
                     limits = c(X_LO, X_HI), expand = c(0, 0)) +
  labs(
    title = "Per-gene change in editing ratio, by gene set",
    # Shortened to two lines. The selection caveat stays: without it the
    # subsets' 100% negativity reads as independent evidence of directional
    # consistency when it is partly a consequence of how they were chosen.
    subtitle = paste(strwrap(paste(
      "One value per gene: mean editing ratio in knockdown minus scrambled control.",
      "Box = IQR, tick = median, diamond = mean, whiskers = 5th-95th; densities are width-scaled, not n-scaled.",
      "The lower three sets are subsets of the first, selected on significance and so on effect size,",
      "which partly builds in their stronger negative shift."),
      width = 158), collapse = "\n"),
    x = "Change in editing ratio (knockdown - scrambled control)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14.5, colour = INK),
    plot.subtitle = element_text(size = 9.1, colour = MUTED, margin = margin(b = 12)),
    strip.text = element_text(face = "bold", size = 11.5, colour = INK, hjust = 0,
                              margin = margin(b = 5, t = 3)),
    axis.title.x = element_text(size = 10, colour = MUTED, margin = margin(t = 9)),
    axis.text.y = element_text(size = 10.2, colour = INK, hjust = 1),
    axis.text.x = element_text(size = 9.2, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = GRIDS, linewidth = 0.35),
    panel.spacing.y = unit(1.5, "lines"),
    plot.background = element_rect(fill = SURF, colour = NA),
    plot.margin = margin(13, 15, 11, 12))

ggsave("delta_distributions_endo2.png", p, width = 11.6, height = 7.4, dpi = 300, bg = SURF)
cat("Written: delta_distributions_endo2.png\n")
print(stat[, .(comp, set, n, pneg = round(pneg, 1), mean = round(mu, 4),
               median = round(med, 4), IQR = round(q3 - q1, 4))])
