# Figure 7. Where editing sites fall within genes, control against knockdown.
#
# Every editing site was annotated with the Ensembl Variant Effect Predictor,
# which reports what part of a transcript the site sits in and what effect a
# base change there would have. Examples are the three prime untranslated
# region, an intron, or a change that alters an amino acid.
#
# The left panel is the scrambled control, the right panel the two knockdowns
# combined. The comparison asks whether knocking down ADAR changes the kind of
# place editing happens, not just how much of it there is.
#
# WHY THIS IS A BAR CHART RATHER THAN A PIE CHART
#
# There are nine categories, and seven of them account for under two per cent of
# sites each. A pie chart cannot show a nine way split at those proportions, so
# the earlier version showed the largest three and collapsed the rest into
# "Other". That hid precisely the categories a reader would want to compare
# between control and knockdown, such as sites that change an amino acid.
#
# WHY THE AXIS IS ON A LOG SCALE
#
# Three prime untranslated regions and introns together hold about 98 per cent
# of all sites. The remaining seven categories range from 2 sites up to 126. On
# an ordinary axis those seven bars would be too short to see. A log scale gives
# each category a visible bar.
#
# The cost is that bar lengths can no longer be compared by eye as ratios, since
# equal spacing on a log axis means equal multiplication rather than equal
# addition. The exact count and percentage are printed at the end of every bar
# so the real numbers are always available.
#
# WHY BOTH PANELS USE THE SAME CATEGORY ORDER
#
# Categories are ordered once, by how common they are across control and
# knockdown combined, and that order is then fixed in both panels. If each panel
# sorted itself, a category could sit in a different row on each side and the
# two panels could not be compared by position.
#
# WHY EACH PANEL IS A SINGLE COLOUR
#
# Colour distinguishes the two conditions and nothing else. Giving each
# consequence category its own colour would add a nine entry legend that
# duplicates the row labels already on the axis.
#
# Input:  vep_site_dedup.txt (knockdown), vep_site_dedup_scr.txt (control)
# Output: consequence_bar.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

OUT <- "/rds/general/user/sj1825/home/endothelial2/output"
vep_ko  <- fread(file.path(OUT, "vep_site_dedup.txt"))
vep_scr <- fread(file.path(OUT, "vep_site_dedup_scr.txt"))

# VEP writes category names with underscores, for example three_prime_UTR. This
# swaps them for spaces so the axis labels read as ordinary English.
shorten <- function(x) gsub("_", " ", x)

# Fix the category order once, using both conditions pooled, so the two panels
# stay aligned row for row. The order is reversed because the panels are turned
# on their side below, and without reversing, the most common category would end
# up at the bottom.
pooled <- rbind(vep_ko[, .(primary)], vep_scr[, .(primary)])[, .N, by = primary][order(-N)]
term_order <- rev(pooled$primary)

# Counts one panel's sites by category. Both panels go through this same
# function so they are built identically.
make_panel_data <- function(vep_gb) {
  d <- vep_gb[, .N, by = primary]
  # Join onto the full fixed category list. A category present in one condition
  # but absent from the other would otherwise disappear from that panel, and the
  # two panels would no longer line up. Missing categories become an explicit
  # count of zero instead.
  d <- merge(data.table(primary = term_order), d, by = "primary", all.x = TRUE)
  d[is.na(N), N := 0L]
  d[, primary := factor(shorten(primary), levels = shorten(term_order))]
  d[, pct := N / sum(N)]
  d
}

d_scr <- make_panel_data(vep_scr)
d_ko  <- make_panel_data(vep_ko)

SCR_COLOR <- "#8a897f"
ACCENT <- "#2a78d6"
TEXT_MUTED <- "#898781"

bar_panel <- function(d, color, title, subtitle) {
  ggplot(d, aes(x = primary, y = N)) +
    geom_col(fill = color, width = 0.68) +
    # Print the count and percentage at the end of each bar. This matters more
    # than usual here, because the log scale makes bar lengths hard to compare
    # by eye. A zero count has no bar at all, so it is labelled explicitly to
    # show the category was looked for and not found.
    geom_text(aes(label = ifelse(N == 0, "0", sprintf("%s (%.1f%%)", format(N, big.mark = ","), 100 * pct))),
              hjust = -0.08, size = 3.3, color = "#0b0b0b") +
    coord_flip(clip = "off") +
    scale_y_log10(expand = expansion(mult = c(0, 0.55)), labels = scales::label_number(big.mark = ",")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Sites (log scale)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "#52514e", size = 9, hjust = 0.5, margin = margin(b = 10)),
      axis.text.y = element_text(size = 10.5, color = "#0b0b0b"),
      axis.text.x = element_text(color = TEXT_MUTED, size = 8.5),
      axis.title.x = element_text(color = "#52514e", size = 9, margin = margin(t = 6)),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "#e1e0d9", linewidth = 0.4),
      plot.margin = margin(10, 40, 10, 10)
    )
}

p_scr <- bar_panel(d_scr, SCR_COLOR, "scr (baseline)", sprintf("%s gene-attached sites", format(nrow(vep_scr), big.mark = ",")))
p_ko  <- bar_panel(d_ko, ACCENT, "siADAR1 + siADAR2 knockdown", sprintf("%s gene-attached sites", format(nrow(vep_ko), big.mark = ",")))

combined <- (p_scr | p_ko) +
  plot_annotation(
    title = "Variant consequence types (all 9 categories): scr baseline vs ADAR1/ADAR2 knockdown",
    subtitle = "Primary consequence term per site, log-scaled x-axis. Same category order both panels.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle = element_text(color = "#52514e", size = 9.5, hjust = 0.5, margin = margin(b = 6))
    )
  ) &
  theme(plot.background = element_rect(fill = "#fcfcfb", color = NA))

ggsave(file.path(OUT, "consequence_bar.png"), combined, width = 12.5, height = 5.8, dpi = 300, bg = "#fcfcfb")
cat("Written:", file.path(OUT, "consequence_bar.png"), "\n")
cat("\n-- scr --\n"); print(d_scr[, .(primary, N, pct = sprintf("%.2f%%", 100*pct))])
cat("\n-- knockdown --\n"); print(d_ko[, .(primary, N, pct = sprintf("%.2f%%", 100*pct))])
