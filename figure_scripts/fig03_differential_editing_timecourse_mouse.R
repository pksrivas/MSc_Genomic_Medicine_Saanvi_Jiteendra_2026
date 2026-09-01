# Figure 3. Change in RNA editing across the dehydration time course, mouse.
#
# Mice were dehydrated and kidney endothelial cells sampled at five timepoints,
# from 0 hours (the control) to 48 hours. Cells come from three compartments of
# the kidney vasculature: MEC, CEC and GEC. Each compartment has its own 0 hour
# control library.
#
# This figure plots the CHANGE in editing rather than the level of it. For every
# library, editing is expressed as the difference from that compartment's own 0
# hour control, so MEC is compared against MEC1, CEC against CEC1 and GEC
# against GEC1. A point above the dashed zero line means more editing than at
# the start, below it means less.
#
# Comparing each compartment against its own control rather than against a
# shared baseline matters because the three compartments do not start at the
# same editing level. A shared baseline would show differences between
# compartments that were already there before dehydration began.
#
# WHAT THIS FIGURE IS NOT
#
# These are differences between two averages. They have not been put through a
# significance test, and no p values are attached to them. The figure describes
# what the data looks like, and no claim of statistical significance should be
# read into a point sitting away from zero.
#
# A per timepoint significance test was deliberately not built, because the
# design cannot support one. Each combination of compartment and timepoint is a
# single sequencing library. A test comparing one timepoint against control
# would have no repeated libraries within that timepoint to measure variation
# against, so any p value it produced would reflect variation between cell
# clusters within one animal rather than between animals. That is the same
# pseudoreplication problem the rest of this project takes care to avoid. A
# genuine per timepoint test would need building explicitly, and is flagged here
# rather than quietly approximated.
#
# THE TWO MISSING CEC TIMEPOINTS ARE MISSING FOR DIFFERENT REASONS
#
# CEC contributes no data at 24 or 36 hours, but not for the same reason, and
# the two should not be described together.
#
# At 24 hours the library exists. CEC3 was sequenced and is deposited in
# E-MTAB-8145. It is absent here because none of its cells survived the
# published cluster assignment, so there were no cells to aggregate into
# pseudobulk units.
#
# At 36 hours there is no library at all. No CEC sample was generated for that
# timepoint, so nothing was ever sequenced.
#
# The first is a loss during processing, the second a gap in the experiment.
#
# HOW THE ERROR BARS ARE CALCULATED
#
# Each point is a difference between two averages, and each of those averages
# carries its own uncertainty. Because the timepoint library and the control
# library are made of different cells, the two uncertainties are independent and
# combine as the square root of the sum of their squares. The error bars show
# one standard error either side.
#
# Input:  all_ec_editing.txt, the per site per pseudobulk unit count table
#         produced by the SPRINT detection stage
# Output: differential_editing_over_time_mouse.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

EPH <- "/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full"
OUT <- "/rds/general/user/sj1825/home/msc_prj/test_data_mouse"

# Read the count table. Each row is one editing site in one pseudobulk unit,
# where a unit is one cell cluster within one sequencing library.
raw <- fread(file.path(EPH, "all_ec_editing.txt"))

# Sample names look like MEC1_cluster_3, so stripping the cluster suffix gives
# the library, and the first three characters of that give the compartment.
raw[, library := sub("_cluster_.*", "", sample)]
raw[, compartment := substr(library, 1, 3)]

# Map each library to its dehydration timepoint. The numbering is consistent
# across compartments, so library 1 is always the 0 hour control and library 5
# is always 48 hours. CEC is absent at 24 and 36 hours; see the header for why
# those two gaps have different causes.
raw[, timepoint := fcase(
  library %in% c("MEC1", "CEC1", "GEC1"), "0h",
  library %in% c("MEC2", "CEC2", "GEC2"), "12h",
  library %in% c("MEC3", "GEC3"),         "24h",
  library %in% c("MEC4", "GEC4"),         "36h",
  library %in% c("MEC5", "CEC5", "GEC5"), "48h"
)]
# Stop immediately if any library failed to match the list above. Without this
# an unrecognised library would silently become a missing timepoint and drop out
# of the figure without anyone noticing.
stopifnot(!anyNA(raw$timepoint))
raw[, timepoint := factor(timepoint, levels = c("0h", "12h", "24h", "36h", "48h"))]

# Keep only sites covered by at least 10 reads. Editing ratios from very low
# coverage are unstable, since a single read changes them a great deal.
pass <- raw[total >= 10]

# Collapse to one editing ratio per pseudobulk unit. Reads are summed across all
# that unit's sites first and the ratio taken afterwards, so the result is the
# proportion of reads that were edited rather than an average of per site
# proportions. Averaging the proportions would give a site with 10 reads the
# same weight as one with 1,000.
unit_ratio <- pass[, .(edited = sum(edited), total = sum(total)), by = .(sample, library, compartment, timepoint)]
unit_ratio[, ratio := edited / total]

# Average those unit level ratios up to one value per library, and record the
# standard error, which measures how much the units within a library disagree.
lib_summary <- unit_ratio[, .(
  n_units = .N,
  mean_ratio = mean(ratio),
  se_ratio = sd(ratio) / sqrt(.N)
), by = .(library, compartment, timepoint)]

# ---------------------------------------------------------------------------
# Express every timepoint as a change from its own compartment's control
# ---------------------------------------------------------------------------

# Pull out the 0 hour library for each compartment, then join it back onto the
# remaining timepoints by compartment so every row carries its own baseline.
ctrl <- lib_summary[timepoint == "0h", .(compartment, ctrl_mean = mean_ratio, ctrl_se = se_ratio, ctrl_n = n_units)]
delta_dt <- merge(lib_summary[timepoint != "0h"], ctrl, by = "compartment")

# The change itself, and its uncertainty. The two libraries contain different
# cells, so their errors are independent and combine in quadrature.
delta_dt[, delta := mean_ratio - ctrl_mean]
delta_dt[, se_delta := sqrt(se_ratio^2 + ctrl_se^2)]
setorder(delta_dt, compartment, timepoint)
print(delta_dt[, .(compartment, timepoint, mean_ratio, ctrl_mean, delta, se_delta, n_units, ctrl_n)])

all_levels <- c("12h", "24h", "36h", "48h")
delta_dt[, timepoint := factor(as.character(timepoint), levels = all_levels)]

PAL <- c(MEC = "#2a78d6", CEC = "#eb6834", GEC = "#1baf7a")

# CEC has no data at 24 or 36 hours. Drawing a plain line through its points
# would join 12 hours straight to 48 hours and look identical to a line through
# measured data, implying the trend between them was observed when it was not.
#
# This block builds the connecting lines by hand instead, one segment at a time,
# and marks any segment that skips a timepoint. Those are drawn dashed further
# down, so the gap in the sampling is visible in the figure.
seg_dt <- rbindlist(lapply(split(delta_dt, delta_dt$compartment), function(d) {
  d <- d[order(timepoint)]
  if (nrow(d) < 2) return(NULL)
  rbindlist(lapply(seq_len(nrow(d) - 1), function(i) {
    lvl_from <- which(all_levels == as.character(d$timepoint[i]))
    lvl_to   <- which(all_levels == as.character(d$timepoint[i + 1]))
    data.table(compartment = d$compartment[1],
               x = d$timepoint[i], xend = d$timepoint[i + 1],
               y = d$delta[i], yend = d$delta[i + 1],
               is_gap = (lvl_to - lvl_from) > 1)
  }))
}))

p <- ggplot(delta_dt, aes(x = timepoint, y = delta, color = compartment, group = compartment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#898781", linewidth = 0.5) +
  geom_segment(data = seg_dt, aes(x = x, xend = xend, y = y, yend = yend, linetype = is_gap),
               linewidth = 0.9, inherit.aes = TRUE) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  geom_errorbar(aes(ymin = delta - se_delta, ymax = delta + se_delta), width = 0.08, linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(size = n_units)) +
  scale_color_manual(values = PAL, name = "Compartment") +
  scale_size_continuous(range = c(2, 4.5), name = "Pseudobulk\nunits (clusters)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Differential RNA editing across the dehydration time-course, by compartment",
    subtitle = paste(strwrap(
      "Error bars = propagated standard error (SE) from both libraries' own pseudobulk-unit spread. Dashed segment (CEC) = no library at 24h/36h",
      width = 130), collapse = "\n"),
    x = "Dehydration timepoint", y = "Delta editing ratio (timepoint - own 0h control)"
  ) +
  theme_minimal(base_size = 12.5) +
  theme(
    plot.title = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 9.3, color = "#52514e", margin = margin(b = 12)),
    axis.text = element_text(color = "#52514e"),
    axis.title = element_text(color = "#52514e", size = 10.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#e1e0d9", linewidth = 0.4),
    legend.position = "right",
    plot.background = element_rect(fill = "#fcfcfb", color = NA),
    plot.margin = margin(12, 16, 12, 12)
  )

ggsave(file.path(OUT, "differential_editing_over_time_mouse.png"), p, width = 9.8, height = 5.8, dpi = 300, bg = "#fcfcfb")
cat("Written:", file.path(OUT, "differential_editing_over_time_mouse.png"), "\n")
