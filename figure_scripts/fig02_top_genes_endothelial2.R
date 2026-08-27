# Figure 2. Genes with the most differentially edited sites, endothelial2.
#
# One horizontal bar per gene, showing how many of that gene's editing sites the
# GLMM called significant after knockdown. The left panel is ADAR1 knockdown
# against scrambled control, the right panel ADAR2 knockdown against the same
# control.
#
# WHAT THE RANKING IS BASED ON
#
# Genes are ranked by the number of sites that came out SIGNIFICANT, not by the
# number of sites tested. The distinction matters. A long, heavily sequenced
# gene can carry hundreds of tested sites and show no knockdown effect at any of
# them. Ranking by tested sites would put such a gene at the top of a figure
# that sits directly under a discussion of knockdown effects, which invites the
# reader to draw a conclusion the data does not support.
#
# The two comparisons are also ranked separately rather than pooled. Pooling
# them produces numbers that belong to neither comparison. CTSB, for example,
# has 136 significant sites under siADAR1 and 138 under siADAR2, but a pooled
# count of 156, which describes neither.
#
# WHY THERE IS NO THIRD PANEL FOR THE CONTROL
#
# The scrambled control is the baseline both knockdowns are compared against. It
# is not compared against itself, so the phrase "differentially edited sites"
# has no meaning for it and there is no third bar chart to draw. The control is
# shown in its own right in the consequence and raw editing ratio figures.
#
# WHY THE RIGHT PANEL HAS ONLY FOUR BARS
#
# ADAR2 knockdown produced only four genes with any significant site at all, and
# each of those has exactly one. The panel shows four bars because there are
# four genes, not because it has been truncated. Padding it out to ten to match
# the left panel would hide how little that comparison found.
#
# Input:  gene_table_siADAR1.txt, gene_table_siADAR2.txt
# Output: top_genes_bar.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUT <- "/rds/general/user/sj1825/home/endothelial2/output"
g1 <- fread(file.path(OUT, "gene_table_siADAR1.txt"))
g2 <- fread(file.path(OUT, "gene_table_siADAR2.txt"))

ACCENT_1 <- "#2a78d6"
ACCENT_2 <- "#eb6834"
TEXT_MUTED <- "#898781"

# The axis on these panels counts sites, so it should only ever be labelled with
# whole numbers. R's default axis break function does not know that. On the
# siADAR2 panel, where the largest bar is 1, it would label the axis 0, 0.2,
# 0.4, 0.6, 0.8, 1, which reads as though a gene could have a fifth of a site.
# This wrapper keeps only the whole number breaks.
integer_breaks <- function(n = 5) function(x) { b <- pretty(x, n); unique(b[b == floor(b)]) }

# Draws one panel. The same function is used for both comparisons so they are
# guaranteed to be built identically and can be compared by eye.
bar_panel <- function(g, color, comparison_label, n_max = 10L) {
  # Keep genes with at least one significant site, put the largest first, and
  # take the top n_max. If fewer than n_max genes qualify, all of them are kept
  # and the panel simply has fewer bars, which is what happens for siADAR2.
  d <- g[n_GLMM_sig >= 1][order(-n_GLMM_sig)]
  d <- head(d[, .(SYMBOL, n_GLMM_sig)], n_max)
  # Fixing the gene order as a factor stops ggplot resorting them alphabetically.
  # The order is reversed because the panel is flipped on its side below, and
  # flipping would otherwise put the largest bar at the bottom.
  d[, SYMBOL := factor(SYMBOL, levels = rev(SYMBOL))]

  ggplot(d, aes(x = SYMBOL, y = n_GLMM_sig)) +
    geom_col(fill = color, width = 0.68) +
    geom_text(aes(label = n_GLMM_sig), hjust = -0.3, size = 4.0, fontface = "bold", color = "#0b0b0b") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)), breaks = integer_breaks()) +
    labs(title = comparison_label,
         subtitle = sprintf("%d genes with >=1 GLMM-significant site (of %s tested)",
                             sum(g$n_GLMM_sig >= 1), format(nrow(g), big.mark = ",")),
         x = NULL, y = "GLMM-significant sites (FDR<0.05)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "#52514e", size = 9, hjust = 0.5, margin = margin(b = 10)),
      axis.text.y = element_text(size = 11, face = "italic", color = "#0b0b0b"),
      axis.text.x = element_text(color = TEXT_MUTED),
      axis.title.x = element_text(color = "#52514e", size = 9.5, margin = margin(t = 6)),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#e1e0d9", linewidth = 0.4),
      plot.margin = margin(10, 22, 10, 10)
    )
}

p1 <- bar_panel(g1, ACCENT_1, "siADAR1", n_max = 10L)
p2 <- bar_panel(g2, ACCENT_2, "siADAR2", n_max = 10L)

combined <- (p1 | p2) +
  plot_annotation(
    # The overall subtitle has been removed from the figure and its content
    # moved into the written caption. It explained two things a reader cannot
    # work out from the bars alone: that the ranking uses significant sites
    # rather than tested sites, and that siADAR2's four bars are the complete
    # set rather than a shortened list. Both still need saying, just not here.
    #
    # The per panel subtitles stay, because they carry the denominators, 192 of
    # 725 genes for siADAR1 and 4 of 467 for siADAR2. Those are measurements
    # rather than explanation, and the panels are misleading without them.
    title = "Top genes by differentially edited (GLMM-significant) site count",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0, margin = margin(b = 12))
    )
  ) &
  theme(plot.background = element_rect(fill = "#fcfcfb", color = NA))

ggsave(file.path(OUT, "top_genes_bar.png"), combined, width = 11, height = 5.3, dpi = 300, bg = "#fcfcfb")
cat("Written:", file.path(OUT, "top_genes_bar.png"), "\n")

cat("\n-- siADAR1 top 10 by n_GLMM_sig --\n")
print(head(g1[n_GLMM_sig >= 1][order(-n_GLMM_sig), .(SYMBOL, n_GLMM_sig, n_sites)], 10))
cat("\n-- siADAR2, all genes with n_GLMM_sig >= 1 --\n")
print(g2[n_GLMM_sig >= 1][order(-n_GLMM_sig), .(SYMBOL, n_GLMM_sig, n_sites)])
