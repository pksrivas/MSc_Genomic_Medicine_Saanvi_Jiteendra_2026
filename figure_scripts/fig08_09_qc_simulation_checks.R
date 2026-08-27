# Figures 8 and 9. Quality checks on the simulated data.
#
# Figure 6 uses simulated data to compare the three statistical tests. That
# argument only holds if the simulator does what it claims. These two figures
# are the check on it, and this one script produces both.
#
#   qc_parameter_recovery.png  Figure 8. Did the simulator produce data with the
#                              properties it was asked for?
#
#   qc_dispersion.png          Figure 9. Is the simulated data as messy as real
#                              data, or is it an easier problem than reality?
#
# The second question is the more important of the two. A simulator that
# produces unrealistically clean data would make every test look better behaved
# than it is, and conclusions drawn from it would not transfer.
#
# FIGURE 8, IN FOUR PANELS
#
# Each panel asks for one property to be recovered. The dashed line in every
# panel is where the points should land if the simulator is correct, so the
# check is simply whether the points sit on the line.
#
#   A  Effect sizes. Three differences were planted, of 0.05, 0.10 and 0.20.
#      This plots what was asked for against what came out.
#   B  Baseline editing rate, which was set to 0.10.
#   C  Mean read coverage, which was set to 30.
#   D  Between-sample variation. This is sigma, the same quantity varied along
#      the horizontal axis of Figure 6, checked against the value requested.
#
# Colour is the sigma setting and the filled or open symbol is the sample size,
# so each point can be traced back to the configuration that produced it.
#
# FIGURE 9, AND WHAT DISPERSION MEANS
#
# Dispersion measures how much more variable the counts are than a simple
# binomial model expects. A value of 1 means the data is exactly as variable as
# that model predicts. Above 1 means samples disagree more than the model allows
# for, which is the situation the random effect in the GLMM exists to handle.
#
# Simulated datasets at three sigma settings are shown in blue, and the four
# real datasets in orange. If the orange distributions sat far above the blue
# ones, the simulation would be an easier problem than reality.
#
# TWO THINGS ABOUT HOW FIGURE 9 IS DRAWN
#
# The vertical axis is logarithmic. Dispersion spans an enormous range, and on
# an ordinary axis every distribution would be squashed against the bottom.
#
# A log axis cannot show a value of zero, and it cannot show the extreme tails
# without making everything else unreadable, so the drawn shapes exclude zeros
# and the axis is clipped. The medians plotted on top are deliberately
# calculated from the complete distribution before any of that, so they match
# the reported table exactly rather than the trimmed shape. The proportion of
# values excluded is printed to the console and belongs in the figure caption.
#
# Input:  qc_reditr_simulations.txt for Figure 8, plus the four observed count
#         tables and their metadata for Figure 9
# Output: qc_parameter_recovery.png, qc_dispersion.png

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(reditR)
})

H <- "/rds/general/user/sj1825/home"
setwd(file.path(H, "msc_prj"))

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; PURPLE <- "#7a4ea8"
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; MUTED <- "#52514e"; GRIDS <- "#e1e0d9"

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 12.5, colour = INK,
                                margin = margin(b = 10)),
    axis.title   = element_text(size = 10, colour = MUTED),
    axis.text    = element_text(size = 9,  colour = MUTED),
    legend.position = "top", legend.justification = "left",
    legend.text  = element_text(size = 10, colour = INK),
    legend.title = element_text(size = 10, colour = MUTED),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = GRIDS, linewidth = 0.35),
    plot.background  = element_rect(fill = SURF, colour = NA),
    plot.margin = margin(10, 12, 10, 10))

# ===========================================================================
# Figure 8: did the simulator produce what it was asked for?
# ===========================================================================

# The first six rows of the QC table are the parameter recovery section, one row
# per simulator configuration. Later rows hold the dispersion comparison, which
# this figure does not use.
rec <- fread("qc_reditr_simulations.txt", nrows = 6)
rec[, `:=`(sigma_f = factor(sigma_set, levels = c(0, 0.25, 1),
                            labels = c("0", "0.25", "1.0")),
           n_f     = factor(npc_set,  levels = c(3, 6),
                            labels = c("n = 3", "n = 6")))]

sig_cols <- c("0" = BLUE, "0.25" = ORANGE, "1.0" = PURPLE)

# Panel A. Three effect sizes were planted in every configuration, so each
# configuration contributes three points. Reshaping to one row per planted
# effect lets all three be plotted against what was requested.
eff <- melt(rec, id.vars = c("sigma_f", "n_f"),
            measure.vars = c("delta_05", "delta_10", "delta_20"),
            variable.name = "which", value.name = "observed")
eff[, nominal := c(delta_05 = 0.05, delta_10 = 0.10, delta_20 = 0.20)[as.character(which)]]

pA <- ggplot(eff, aes(nominal, observed, colour = sigma_f, shape = n_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = MUTED, linewidth = 0.5) +
  geom_point(size = 2.6, stroke = 0.9, fill = NA) +
  scale_colour_manual(values = sig_cols, name = "sigma") +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  scale_x_continuous(breaks = c(0.05, 0.10, 0.20), limits = c(0.02, 0.23)) +
  scale_y_continuous(limits = c(0.02, 0.23)) +
  labs(title = "A  Planted effect recovered", x = "Nominal effect size", y = "Observed median") +
  base_theme

# B: baseline editing rate
pB <- ggplot(rec, aes(sigma_set, baseline_obs, colour = sigma_f, shape = n_f)) +
  geom_hline(yintercept = 0.10, linetype = "22", colour = MUTED, linewidth = 0.5) +
  geom_point(size = 2.6, stroke = 0.9, fill = NA) +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_shape_manual(values = c(16, 1), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 1.0)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(title = "B  Baseline editing rate", x = "sigma set", y = "Observed baseline") +
  base_theme

# C: mean coverage
pC <- ggplot(rec, aes(sigma_set, coverage_obs, colour = sigma_f, shape = n_f)) +
  geom_hline(yintercept = 30, linetype = "22", colour = MUTED, linewidth = 0.5) +
  geom_point(size = 2.6, stroke = 0.9, fill = NA) +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_shape_manual(values = c(16, 1), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 1.0)) +
  scale_y_continuous(limits = c(29, 31)) +
  labs(title = "C  Mean coverage", x = "sigma set", y = "Observed coverage") +
  base_theme

# D: logit-scale SD
pD <- ggplot(rec, aes(sigma_set, sd_logit_obs, colour = sigma_f, shape = n_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = MUTED, linewidth = 0.5) +
  geom_point(size = 2.6, stroke = 0.9, fill = NA) +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_shape_manual(values = c(16, 1), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 1.0), limits = c(-0.05, 1.05)) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "D  Between-sample SD, logit scale", x = "sigma set", y = "Observed SD") +
  base_theme

fig1 <- (pA | pB) / (pC | pD) +
  plot_annotation(theme = theme(plot.background = element_rect(fill = SURF, colour = NA)))
ggsave("qc_parameter_recovery.png", fig1, width = 10.5, height = 8.2, dpi = 300, bg = SURF)
cat("Written: qc_parameter_recovery.png\n")

# ===========================================================================
# Figure 9: is the simulated data as variable as the real data?
# ===========================================================================

# Dispersion for one site in one condition, measured as Pearson's chi-square
# statistic divided by its degrees of freedom.
#
# In words: for each sample, compare the number of edited reads actually seen
# against the number a binomial model would predict from the pooled editing
# rate, scale that difference by how much the model says it should vary, and
# average across samples. The result is 1 when the data behaves exactly as the
# model expects, and larger when samples disagree more than it allows.
#
# Two situations return a missing value rather than a number. A site seen in
# only one sample gives nothing to compare against. A site where every read is
# edited, or none are, has no variance for the model to predict.
dispersion <- function(d) {
  d <- d[total > 0]
  d[, { n <- .N
    if (n < 2L) NA_real_ else { ph <- sum(edited)/sum(total)
      if (ph <= 0 || ph >= 1) NA_real_ else
        sum((edited - total*ph)^2/(total*ph*(1-ph)))/(n-1) } }, by = .(site, condition)]$V1
}

# Simulate at three sigma settings and keep only the sites with no planted
# effect. Sites carrying a real difference between conditions would show extra
# variation for that reason alone, which is not what is being measured here.
#
# The seed is fixed so the figure is reproducible.
set.seed(20260812)
sim <- rbindlist(lapply(c(0, 0.25, 1.0), function(sg) {
  s <- simulate_editing_data(n_null = 800L, n_effects = c(`0.10` = 100),
                             n_per_condition = 6L, sample_re_sd = sg, seed = 7L)
  e <- merge(as.data.table(s$editing), as.data.table(s$metadata), by = "sample")
  e <- merge(e, as.data.table(s$truth), by = "site")[true_effect == 0]
  data.table(dataset = sprintf("Simulated\nsigma = %.2f", sg), class = "Simulated",
             phi = dispersion(e))
}))

# The four observed datasets, each reduced to the same five columns so they can
# go through the identical dispersion calculation.
real <- list()
dia <- merge(fread(file.path(H, "diabetes_output/filtered_sites_clustered_t.txt")),
             fread(file.path(H, "diabetes_output/sample_metadata.txt")), by = "sample")
real[["Cardiomyocyte"]] <- dia[, .(site, sample, edited, total, condition)]
# Read from $HOME rather than ephemeral storage. The ephemeral copy of this
# table is liable to be cleared without warning, which already happened to
# another mouse input file used elsewhere in this folder.
mo <- fread(file.path(H, "msc_prj/test_data_mouse/sprint_output_full/all_ec_clustered_with_condition.txt"))
real[["Mouse\npseudobulk"]] <- mo[, .(site, sample, edited, total, condition)]
en <- merge(fread(file.path(H, "endothelial2/output/filtered_sites_clustered.txt")),
            fread(file.path(H, "endothelial2/output/sample_metadata.txt")), by = "sample")
en <- en[, .(site, sample, edited, total, condition)]
real[["Endothelial\nall sites"]] <- en
# The endothelial data appears twice, once complete and once restricted to
# sites covered in every sample. Sites missing from some samples are measured
# on fewer observations, so including them could make the dataset look more or
# less variable than it is. Showing both separates that effect from the
# underlying comparison.
cplt <- en[total > 0, uniqueN(sample), by = site][V1 == uniqueN(en$sample), site]
real[["Endothelial\ncomplete sites"]] <- en[site %in% cplt]

obs <- rbindlist(lapply(names(real), function(nm)
  data.table(dataset = nm, class = "Observed", phi = dispersion(real[[nm]]))))

all_phi <- rbind(sim, obs)[!is.na(phi)]
all_phi[, dataset := factor(dataset, levels = c(unique(sim$dataset), names(real)))]

# Order matters in the next three lines. The medians are calculated first, from
# the complete distribution, so they match the numbers reported in the text.
# Only afterwards are zero values removed for drawing, because a log axis cannot
# place them. Calculating the medians after that removal would shift them, which
# is a mistake worth guarding against explicitly.
#
# The proportion of values dropped, and the proportion falling outside the
# visible range, are recorded here and printed at the end. Both belong in the
# figure caption so a reader knows what the drawn shapes leave out.
meds <- all_phi[, .(m = median(phi)), by = .(dataset, class)]
zero_frac <- all_phi[, .(pct_zero    = round(100 * mean(phi == 0), 1),
                         pct_clipped = round(100 * mean(phi < 0.01 | phi > 30), 1)),
                     by = dataset]
dt <- all_phi[phi > 0]

fig2 <- ggplot(dt, aes(dataset, phi, fill = class, colour = class)) +
  geom_hline(yintercept = 1, linetype = "22", colour = MUTED, linewidth = 0.5) +
  geom_violin(alpha = 0.30, linewidth = 0.4, scale = "width", trim = TRUE) +
  geom_point(data = meds, aes(dataset, m), size = 2.6, stroke = 0, show.legend = FALSE) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10),
                labels = c("0.01", "0.1", "1", "10")) +
  coord_cartesian(ylim = c(0.01, 30)) +
  scale_fill_manual(values = c(Simulated = BLUE, Observed = ORANGE), name = NULL) +
  scale_colour_manual(values = c(Simulated = BLUE, Observed = ORANGE), name = NULL) +
  labs(title = "Per-site dispersion, simulated against observed data",
       x = NULL, y = "Dispersion (Pearson chi-square / df), log scale") +
  base_theme +
  theme(axis.text.x = element_text(size = 9, colour = INK, lineheight = 0.95))

ggsave("qc_dispersion.png", fig2, width = 11, height = 6, dpi = 300, bg = SURF)
cat("Written: qc_dispersion.png\n")

print(merge(meds[, .(dataset, class, median_phi = round(m, 2))], zero_frac, by = "dataset")[
  , .(dataset = gsub("\n", " ", dataset), class, median_phi, pct_zero, pct_clipped)])
