# Figure 1. Differential editing testing summary for the endothelial2 dataset.
#
# The experiment knocked down ADAR1 and ADAR2 in endothelial cells and compared
# each against a scrambled control (scr). There are three replicates in every
# group, so each comparison is three samples against three.
#
# The figure has three panels, read left to right:
#
#   Panel 1  How many editing sites each statistical test called significant.
#   Panel 2  How the genes picked out by two of those tests overlap.
#   Panel 3  How many sites were seen at all in each condition, before any
#            statistical testing.
#
# WHY PANEL 1 LABELS EVERY BAR WITH ITS OWN DENOMINATOR
#
# The three tests do not all work on the same sites, so "824 significant sites"
# means nothing until you know out of how many. Fisher's exact test can be run
# at every site that was tested. The GLMM is a more complicated model and does
# not always converge, so it produces a result at fewer sites. The Wilcoxon
# test needs the editing ratios to differ between samples, so it fails at sites
# where every sample gives the same value.
#
# Each bar is therefore labelled "n of evaluable". Written that way, the GLMM's
# 824 calls read as 34.6 per cent of the sites it could actually fit, instead of
# 10.3 per cent of every site in the dataset, which would understate it.
#
# The Wilcoxon bars come out at zero. That is not a near miss. With three
# samples per group the smallest p value the test can possibly produce is 0.077,
# which fails even before correcting for multiple testing.
#
# WHY PANEL 2 USES CIRCLES RATHER THAN A TABLE
#
# The two circles are drawn so their areas are proportional to the number of
# genes in each set, and the overlap is solved numerically to match the number
# of shared genes. A reader can see the relative sizes without reading numbers.
#
# The overlap counts genes that have at least one GLMM-significant site and at
# least one Fisher-significant site. Those need not be the same site, which is
# why the figure says so directly underneath.
#
# Each comparison has its own set of tested genes: 725 for siADAR1 and 467 for
# siADAR2. These are not interchangeable. The figure of 778 that appears in the
# VEP summary is the two sets combined, counting each gene once, so it describes
# the experiment as a whole and is not a denominator for either comparison on
# its own.
#
# WHY PANEL 3 CARRIES A WARNING IN THE CAPTION
#
# Panel 3 counts sites detected per condition with no statistical test involved.
# The three conditions share many sites, so the three bars overlap and do not
# add up to the experiment-wide total. Sequencing depth also differs a lot
# between conditions, which affects how many sites can be detected at all, so
# the depth figures belong in the caption alongside this panel.
#
# WHY COLOUR MEANS TEST AND NOT CONDITION
#
# Colour is used for one thing only, which test produced the bar. Conditions are
# already labelled along the axis, so colouring them as well would spend colour
# on information the reader already has. The three test colours were checked for
# colour vision deficiency and stay distinguishable.
#
# Input:  DRE_siADAR1_v4.txt, DRE_siADAR2_v4.txt, gene_table_siADAR1.txt,
#         gene_table_siADAR2.txt, filtered_sites_clustered.txt,
#         sample_metadata.txt
# Output: diffedit_summary_endo2.png

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggforce)
})
O <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(O)

GLMM_C <- "#2a78d6"; FISH_C <- "#eb6834"; WILC_C <- "#c3c2b7"
NEUTRAL <- "#8a897f"
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; MUTED <- "#52514e"; GRIDS <- "#e1e0d9"; TXT <- "#898781"

d1 <- fread("DRE_siADAR1_v4.txt"); d2 <- fread("DRE_siADAR2_v4.txt")
g1 <- fread("gene_table_siADAR1.txt"); g2 <- fread("gene_table_siADAR2.txt")

base_thm <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, colour = INK),
        plot.subtitle = element_text(colour = MUTED, size = 9, hjust = 0.5, margin = margin(b = 10)),
        axis.text = element_text(colour = TXT), panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = SURF, colour = NA))

# ---------------------------------------------------------------------------
# Panel 1: significant sites per test
# ---------------------------------------------------------------------------

# How many sites each test could evaluate, listed in the order GLMM, Fisher,
# Wilcoxon. A test that returned no p value at a site could not evaluate it, so
# counting the non-missing p values gives the denominator for that test.
#
# Note this uses wilcox_pvalue_all and not wilcox_pvalue. The second column
# holds a follow up test run only at sites the GLMM already called significant,
# so using it here would make the Wilcoxon denominator far too small.
n_eval <- function(d) c(sum(!is.na(d$glmm_pvalue)), sum(!is.na(d$fisher_pvalue)),
                        sum(!is.na(d$wilcox_pvalue_all)))
# How many sites each test called significant. The GLMM and Fisher results are
# already stored as TRUE or FALSE columns. The Wilcoxon count is worked out here
# by applying the same Benjamini-Hochberg correction the other two used, rather
# than being typed in as a zero, so that it would change if the data ever did.
n_sig  <- function(d) c(sum(d$GLMM_sig, na.rm = TRUE), sum(d$Fisher_sig, na.rm = TRUE),
                        sum(p.adjust(d$wilcox_pvalue_all, "BH") < 0.05, na.rm = TRUE))

sites <- rbindlist(list(
  data.table(comp = "siADAR1 vs scr", test = c("GLMM", "Fisher", "Wilcoxon"),
             n = n_sig(d1), evaluable = n_eval(d1)),
  data.table(comp = "siADAR2 vs scr", test = c("GLMM", "Fisher", "Wilcoxon"),
             n = n_sig(d2), evaluable = n_eval(d2))))
sites[, test := factor(test, levels = c("GLMM", "Fisher", "Wilcoxon"))]

p_sites <- ggplot(sites, aes(comp, n, fill = test)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = sprintf("%s\nof %s", format(n, big.mark = ","),
                                format(evaluable, big.mark = ","))),
            position = position_dodge(width = 0.78), vjust = -0.25, lineheight = 0.95,
            size = 3.1, fontface = "bold", colour = INK) +
  scale_fill_manual(values = c(GLMM = GLMM_C, Fisher = FISH_C, Wilcoxon = WILC_C), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
  labs(title = "Significant sites", x = NULL, y = NULL,
       subtitle = "each bar over the sites that test could evaluate") +
  base_thm + theme(legend.position = "top", legend.justification = "center",
                   panel.grid.major.x = element_blank(),
                   axis.text.x = element_text(size = 10.5, colour = INK))

# ---------------------------------------------------------------------------
# Panel 2: gene set overlap for siADAR1, drawn as an area proportional diagram
# ---------------------------------------------------------------------------

# Count the genes in each set. n_g has at least one GLMM-significant site, n_f
# has at least one Fisher-significant site, and n_b has both.
n_g <- sum(g1$n_GLMM_sig >= 1); n_f <- sum(g1$n_Fisher_sig >= 1)
n_b <- sum(g1$n_GLMM_sig >= 1 & g1$n_Fisher_sig >= 1)
# Turn those counts into circle radii. Area equals pi times radius squared, so
# taking the square root of count divided by pi gives a circle whose area is the
# gene count. This is what makes the two circles visually comparable.
r1 <- sqrt(n_g / pi); r2 <- sqrt(n_f / pi)
# Standard geometry for the area shared by two overlapping circles whose centres
# are distance d apart. The two early returns handle the cases where the circles
# do not touch at all, and where the smaller sits entirely inside the larger.
lens <- function(d) { if (d >= r1 + r2) return(0); if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  r1^2 * acos((d^2 + r1^2 - r2^2) / (2*d*r1)) + r2^2 * acos((d^2 + r2^2 - r1^2) / (2*d*r2)) -
  0.5 * sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2)) }
# There is no formula for the centre distance that produces a required overlap,
# so solve for it numerically. uniroot searches between the two extreme cases
# above for the distance at which the shared area equals the shared gene count.
dd <- uniroot(function(d) lens(d) - n_b, c(abs(r1 - r2) + 1e-6, r1 + r2 - 1e-6))$root
# Left and right edges of each circle, used to place the labels. The GLMM count
# is small enough that its label will not fit inside the crescent, so it is
# written outside and joined to the shape with a short line.
x_lb <- -dd/2 - r1; x_rb <- -dd/2 + r1; x_lo <- dd/2 - r2; x_ro <- dd/2 + r2
circ <- data.table(x0 = c(-dd/2, dd/2), y0 = 0, r = c(r1, r2),
                   set = factor(c("GLMM-sig", "Fisher-sig"), levels = c("GLMM-sig", "Fisher-sig")))
lab_in <- data.table(x = c((x_lo + x_rb)/2, (x_rb + x_ro)/2), y = 0,
                     t = c(as.character(n_b), as.character(n_f - n_b)))
callout <- x_lb - 3.0
lab_out <- data.table(x = callout, y = 0, t = as.character(n_g - n_b))
leader <- data.table(x = callout + 0.9, xend = (x_lb + x_lo)/2 - 0.15, y = 0, yend = 0)
set_lab <- data.table(x = c(x_lb - 3.4, x_ro + 1.1), y = r2 + 2.4, h = c(0, 1),
                      set = factor(c("GLMM-sig", "Fisher-sig"), levels = c("GLMM-sig", "Fisher-sig")))

p_genes <- ggplot() +
  geom_circle(data = circ, aes(x0 = x0, y0 = y0, r = r, fill = set), alpha = 0.5,
              colour = TXT, linewidth = 0.5) +
  geom_segment(data = leader, aes(x = x, xend = xend, y = y, yend = yend), colour = TXT, linewidth = 0.45) +
  geom_text(data = lab_in, aes(x, y, label = t), size = 4.1, fontface = "bold", colour = INK) +
  geom_text(data = lab_out, aes(x, y, label = t), size = 4.1, fontface = "bold", colour = GLMM_C) +
  geom_text(data = set_lab, aes(x = x, y = y, label = set, colour = set, hjust = h),
            size = 3.7, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c(`GLMM-sig` = GLMM_C, `Fisher-sig` = FISH_C)) +
  scale_colour_manual(values = c(`GLMM-sig` = GLMM_C, `Fisher-sig` = FISH_C)) +
  coord_equal(xlim = c(x_lb - 4.8, x_ro + 1.9), ylim = c(-r2 - 4.4, r2 + 4.0), clip = "off") +
  annotate("text", x = 0, y = -r2 - 2.9, size = 3.1, colour = MUTED, lineheight = 1.15,
           label = sprintf(paste0("overlap = >=1 significant site under each test, not necessarily the same site\n",
                                  "siADAR2: %d GLMM-sig, %d Fisher-sig of %d genes"),
                           sum(g2$n_GLMM_sig >= 1), sum(g2$n_Fisher_sig >= 1), nrow(g2))) +
  labs(title = "Genes with >=1 significant site",
       subtitle = sprintf("siADAR1, of %s genes tested", nrow(g1))) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, colour = INK),
        plot.subtitle = element_text(colour = MUTED, size = 9, hjust = 0.5, margin = margin(b = 6)),
        legend.position = "none", plot.background = element_rect(fill = SURF, colour = NA),
        plot.margin = margin(6, 10, 10, 10))

# ---------------------------------------------------------------------------
# Panel 3: how many sites were detected in each condition
# ---------------------------------------------------------------------------

# This panel involves no statistical testing at all. It reads the filtered site
# table, attaches each sample's condition, and counts the distinct sites seen in
# each condition together with the total read depth behind them.
x <- merge(fread("filtered_sites_clustered.txt"), fread("sample_metadata.txt"), by = "sample")
det <- x[, .(sites = uniqueN(site), reads = sum(as.numeric(total))), by = condition]
det[, condition := factor(condition, levels = c("scr", "siADAR1", "siADAR2"))]
setorder(det, condition)

p_det <- ggplot(det, aes(condition, sites)) +
  geom_col(fill = NEUTRAL, width = 0.62) +
  geom_text(aes(label = format(sites, big.mark = ",")), vjust = -0.5,
            size = 3.6, fontface = "bold", colour = INK) +
  # Read counts were previously printed inside the bars and have been removed to
  # keep the panel clean. They matter for interpretation, because scr was
  # sequenced more deeply than either knockdown (roughly 0.66 million reads
  # against 0.43 and 0.20 million), and deeper sequencing finds more sites
  # regardless of biology. Since that information is no longer anywhere on the
  # figure, it has to appear in the caption instead.
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Raw site detection by condition",
       subtitle = sprintf("no significance call involved; of %s sites detected experiment-wide, which the conditions overlap rather than partition",
                          format(uniqueN(x$site), big.mark = ",")), x = NULL, y = NULL) +
  base_thm + theme(panel.grid.major.x = element_blank(),
                   axis.text.x = element_text(size = 10.5, colour = INK))

fig <- (p_sites | p_genes | p_det) +
  plot_annotation(
    title = "Differential editing testing summary: endothelial2 (ADAR1 / ADAR2 knockdown vs scrambled control)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0, margin = margin(b = 14)),
                  plot.background = element_rect(fill = SURF, colour = NA))) &
  theme(plot.background = element_rect(fill = SURF, colour = NA))

ggsave("diffedit_summary_endo2.png", fig, width = 15, height = 5.9, dpi = 300, bg = SURF)
cat("Written: diffedit_summary_endo2.png\n")
cat(sprintf("  sites: siADAR1 GLMM %d Fisher %d | siADAR2 GLMM %d Fisher %d\n",
            sum(d1$GLMM_sig, na.rm=TRUE), sum(d1$Fisher_sig, na.rm=TRUE),
            sum(d2$GLMM_sig, na.rm=TRUE), sum(d2$Fisher_sig, na.rm=TRUE)))
cat(sprintf("  genes: siADAR1 %d/%d/%d of %d | siADAR2 %d/%d of %d\n",
            n_g - n_b, n_b, n_f - n_b, nrow(g1),
            sum(g2$n_GLMM_sig >= 1), sum(g2$n_Fisher_sig >= 1), nrow(g2)))
print(det)
