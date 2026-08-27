# Figure 4. RNA editing level by compartment and timepoint, mouse.
#
# The question this figure asks is simple: does the amount of RNA editing go up
# as dehydration continues, and does that happen in all three compartments of
# the kidney vasculature or only some of them?
#
# Figure 3 is the companion to this one. That figure shows the CHANGE from each
# compartment's own control. This one shows the editing level itself, so the
# compartments can be compared against each other rather than only against their
# own starting point.
#
# WHICH LIBRARY BELONGS TO WHICH TIMEPOINT
#
#   0h (control)  MEC1, CEC1, GEC1
#   12h           MEC2, CEC2, GEC2
#   24h           MEC3, GEC3          no CEC library at this timepoint
#   36h           MEC4, GEC4          no CEC library at this timepoint
#   48h           MEC5, CEC5, GEC5
#
# This mapping was taken from the one place in the project where it is written
# down explicitly, rather than inferred from the general description of the
# experiment as covering 12 to 48 hours.
#
# HOW EDITING LEVEL IS MEASURED, AND WHY THE FILTER DIFFERS HERE
#
# For each pseudobulk unit, edited reads and total reads are summed across all
# of that unit's sites, and the editing level is the ratio of the two.
#
# Sites are kept if they have at least 10 reads of coverage. That is the only
# filter applied. Elsewhere in this project a site also has to carry at least 2
# edited reads before it counts as detected, and that extra condition is
# deliberately not used here.
#
# The reason is that requiring 2 edited reads only keeps sites that already show
# editing. Measuring the average editing level across a set of sites chosen for
# showing editing would push the average up by construction, and would do so
# more in samples with more sites, which is exactly the comparison the figure is
# making. Filtering on coverage alone avoids that circularity.
#
# THE LIMITATION THIS FIGURE CANNOT ESCAPE
#
# Every combination of compartment and timepoint is a single sequencing library,
# and each library pools roughly six mice whose individual identities cannot be
# recovered. There is no repetition across animals and none across libraries at
# a given timepoint.
#
# The error bars therefore show how much the cell clusters within one library
# disagree with each other. They do not show how much one mouse differs from
# another, and they cannot be read as though they did. A line that rises across
# the figure is one library's trajectory in each compartment, not a replicated
# dose response. This caveat is printed on the figure itself, not left to the
# caption, because the figure invites the reader to see a trend.
#
# Input:  all_ec_editing.txt, the per site per pseudobulk unit count table
#         produced by the SPRINT detection stage
# Output: editing_over_time_mouse.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

EPH <- "/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full"
OUT <- "/rds/general/user/sj1825/home/msc_prj/test_data_mouse"

# Each row of the count table is one editing site in one pseudobulk unit, where
# a unit is one cell cluster within one sequencing library.
raw <- fread(file.path(EPH, "all_ec_editing.txt"))

# Sample names look like MEC1_cluster_3. Removing the cluster suffix leaves the
# library, and the first three characters of that give the compartment.
raw[, library := sub("_cluster_.*", "", sample)]
raw[, compartment := substr(library, 1, 3)]
raw[, timepoint := fcase(
  library %in% c("MEC1", "CEC1", "GEC1"), "0h",
  library %in% c("MEC2", "CEC2", "GEC2"), "12h",
  library %in% c("MEC3", "GEC3"),         "24h",
  library %in% c("MEC4", "GEC4"),         "36h",
  library %in% c("MEC5", "CEC5", "GEC5"), "48h"
)]
# Fail loudly if a library did not match the list above. Without this check an
# unrecognised library would quietly become a missing timepoint and vanish from
# the figure with no warning.
stopifnot(!anyNA(raw$timepoint))
raw[, timepoint := factor(timepoint, levels = c("0h", "12h", "24h", "36h", "48h"),
                           labels = c("0h (control)", "12h", "24h", "36h", "48h"))]

# Coverage filter only, for the reason given in the header. No minimum number of
# edited reads is imposed, because that would preselect sites that already show
# editing and inflate the very quantity being plotted.
pass <- raw[total >= 10]

# One editing ratio per pseudobulk unit. Reads are summed first and the ratio
# taken afterwards, so a site with 1,000 reads counts for more than one with 10,
# which is the correct weighting when estimating an overall editing level.
unit_ratio <- pass[, .(edited = sum(edited), total = sum(total)), by = .(sample, library, compartment, timepoint)]
unit_ratio[, ratio := edited / total]

# Average up to one value per library. The standard error here measures spread
# across that library's own cell clusters. It is a measure of how uniform the
# cells within a library are, and is not biological replication, since there is
# only one library per compartment and timepoint.
lib_summary <- unit_ratio[, .(
  n_units = .N,
  mean_ratio = mean(ratio),
  se_ratio = sd(ratio) / sqrt(.N)
), by = .(library, compartment, timepoint)]
setorder(lib_summary, compartment, timepoint)
print(lib_summary)

PAL <- c(MEC = "#2a78d6", CEC = "#eb6834", GEC = "#1baf7a")

# CEC was never sampled at 24 or 36 hours. If the points were joined with an
# ordinary line, CEC's 12 hour point would connect straight to its 48 hour point
# and look exactly like a line drawn through measured data, implying the shape
# between them is known.
#
# The connecting lines are therefore built one segment at a time, and any
# segment that jumps over a timepoint is flagged. Flagged segments are drawn
# dashed below, so a reader can see where data is absent.
all_levels <- levels(lib_summary$timepoint)
seg_dt <- rbindlist(lapply(split(lib_summary, lib_summary$compartment), function(d) {
  d <- d[order(timepoint)]
  if (nrow(d) < 2) return(NULL)
  rbindlist(lapply(seq_len(nrow(d) - 1), function(i) {
    lvl_from <- which(all_levels == as.character(d$timepoint[i]))
    lvl_to   <- which(all_levels == as.character(d$timepoint[i + 1]))
    data.table(compartment = d$compartment[1],
               x = d$timepoint[i], xend = d$timepoint[i + 1],
               y = d$mean_ratio[i], yend = d$mean_ratio[i + 1],
               is_gap = (lvl_to - lvl_from) > 1)
  }))
}))

p <- ggplot(lib_summary, aes(x = timepoint, y = mean_ratio, color = compartment, group = compartment)) +
  geom_segment(data = seg_dt, aes(x = x, xend = xend, y = y, yend = yend, linetype = is_gap),
               linewidth = 0.9, inherit.aes = TRUE) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22"), guide = "none") +
  geom_errorbar(aes(ymin = mean_ratio - se_ratio, ymax = mean_ratio + se_ratio), width = 0.08, linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(size = n_units)) +
  scale_color_manual(values = PAL, name = "Compartment") +
  scale_size_continuous(range = c(2, 4.5), name = "Pseudobulk\nunits (clusters)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "RNA editing level across the dehydration time-course, by compartment",
    subtitle = paste(strwrap(
      "Each point = ONE library's mean editing ratio (pooled edited/total across its own tested sites, coverage>=10, no minimum edited-read count). Error bars = spread across that library's own pseudobulk units (clusters) -- NOT biological replication; there is exactly 1 library per (compartment, timepoint) cell, so this is descriptive, not a replicated dose-response. Dashed segment (CEC) = no library at that timepoint (24h/36h); a straight connector, not real data.",
      width = 130), collapse = "\n"),
    x = "Dehydration timepoint", y = "Mean editing ratio (pooled edited/total)"
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

ggsave(file.path(OUT, "editing_over_time_mouse.png"), p, width = 9.5, height = 5.8, dpi = 300, bg = "#fcfcfb")
cat("Written:", file.path(OUT, "editing_over_time_mouse.png"), "\n")
