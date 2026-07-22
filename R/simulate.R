#' Simulate editing data with known ground truth
#'
#' Generates synthetic per-site editing counts under a mixed-effects generative
#' model with known effect sizes. Useful for power analysis and for validating
#' that \code{\link{differential_editing}} recovers planted effects.
#'
#' @param n_null Number of null (no-effect) sites. Default 800.
#' @param n_effects Named numeric vector where names are effect sizes (as
#'   character strings) and values are the number of sites at that effect size.
#'   Default \code{c("0.05" = 100, "0.10" = 100, "0.20" = 100)}.
#' @param n_per_condition Number of samples per condition. Default 6.
#' @param baseline_rate Baseline editing rate in the reference condition.
#'   Default 0.10.
#' @param mean_coverage Mean read coverage per site per sample drawn from
#'   \code{rnbinom(mu = mean_coverage, size = 5)}, floored at 10. Default 30.
#' @param sample_re_sd Standard deviation of the per-sample random effect on
#'   the logit scale. Default 0.3.
#' @param seed Random seed for reproducibility. Default 42.
#'
#' @return A named list with three \code{data.table} elements:
#'   \describe{
#'     \item{\code{editing}}{Per-site, per-sample read counts.}
#'     \item{\code{metadata}}{Sample metadata with \code{sample} and
#'       \code{condition} columns.}
#'     \item{\code{truth}}{Per-site ground truth with \code{site} and
#'       \code{true_effect} columns.}
#'   }
#'
#' @examples
#' sim <- simulate_editing_data(n_null = 50, n_effects = c("0.10" = 20), seed = 1)
#' str(sim)
#'
#' @export
simulate_editing_data <- function(n_null          = 800L,
                                   n_effects       = c("0.05" = 100, "0.10" = 100, "0.20" = 100),
                                   n_per_condition = 6L,
                                   baseline_rate   = 0.10,
                                   mean_coverage   = 30L,
                                   sample_re_sd    = 0.3,
                                   seed            = 42L) {
  set.seed(seed)

  # Sample metadata
  metadata <- data.table(
    sample    = c(paste0("ctrl_", seq_len(n_per_condition)),
                  paste0("case_", seq_len(n_per_condition))),
    condition = factor(
      c(rep("control", n_per_condition), rep("diabetic", n_per_condition)),
      levels = c("control", "diabetic")
    )
  )

  # Site truth table
  null_sites <- data.table(
    site        = paste0("site_null_", seq_len(n_null)),
    true_effect = 0
  )
  effect_sites <- rbindlist(Map(function(eff_name, n) {
    data.table(
      site        = paste0("site_eff", eff_name, "_", seq_len(n)),
      true_effect = as.numeric(eff_name)
    )
  }, names(n_effects), n_effects))

  truth       <- rbind(null_sites, effect_sites)
  all_samples <- metadata$sample
  conditions  <- as.character(metadata$condition)

  # Generate editing counts per site
  rows <- lapply(seq_len(nrow(truth)), function(i) {
    eff     <- truth$true_effect[i]
    site_id <- truth$site[i]

    sample_re  <- rnorm(length(all_samples), 0, sample_re_sd)
    mean_rates <- ifelse(conditions == "diabetic",
                         baseline_rate + eff,
                         baseline_rate)
    mean_rates <- pmax(0.001, pmin(0.999, mean_rates))
    p          <- plogis(qlogis(mean_rates) + sample_re)

    total  <- pmax(10L, rnbinom(length(all_samples), mu = mean_coverage, size = 5))
    edited <- rbinom(length(all_samples), total, p)

    data.table(
      site       = site_id,
      sample     = all_samples,
      edited     = edited,
      total      = total,
      edit_ratio = edited / total
    )
  })

  list(editing = rbindlist(rows), metadata = metadata, truth = truth)
}


#' Validate differential editing results against simulation ground truth
#'
#' Compares the output of \code{\link{differential_editing}} against the
#' ground truth from \code{\link{simulate_editing_data}}. Reports convergence
#' rate (GLMM), false-positive rate on null sites, and power by effect size --
#' independently for whichever test(s) are present in \code{results} (i.e.
#' whichever were requested via \code{differential_editing(test = ...)}).
#' There is no combined verdict to validate, matching
#' \code{differential_editing()}'s independent-tests design: each test is
#' reported on its own.
#'
#' @param results \code{data.table} returned by \code{differential_editing()}.
#' @param truth \code{data.table} with columns \code{site} and
#'   \code{true_effect}, as returned in \code{simulate_editing_data()$truth}.
#' @param fdr_threshold Retained for API consistency; unused -- the
#'   \verb{<Test>_sig} columns in \code{results} already reflect whichever
#'   FDR threshold \code{differential_editing()} was called with.
#'
#' @return An object of class \code{"reditR_validation"} -- a list with
#'   components \code{convergence} (numeric, GLMM only), \code{fpr} (named
#'   list, one element per test present in \code{results}), and \code{power}
#'   (\code{data.table} with one \verb{power_<Test>} column per test present).
#'
#' @examples
#' sim  <- simulate_editing_data(n_null = 30, n_effects = c("0.10" = 10), seed = 1)
#' \dontrun{
#' res  <- differential_editing(tempfile(), tempfile())  # use sim$editing / sim$metadata
#' val  <- validate_against_truth(res, sim$truth)
#' print(val)
#' }
#'
#' @export
validate_against_truth <- function(results, truth, fdr_threshold = 0.05) {
  merged <- merge(results, truth, by = "site", all.x = TRUE)

  convergence <- if ("glmm_pvalue" %in% names(merged))
    mean(!is.na(merged$glmm_pvalue)) else NA_real_

  possible_tests <- c(GLMM = "GLMM_sig", Fisher = "Fisher_sig", Wilcoxon = "Wilcox_sig")
  present_tests  <- possible_tests[possible_tests %in% names(merged)]

  null_sites <- merged[true_effect == 0]
  fpr <- lapply(present_tests, function(cn) {
    if (nrow(null_sites) == 0) return(NA_real_)
    mean(null_sites[[cn]], na.rm = TRUE)
  })
  names(fpr) <- names(present_tests)

  effect_sizes <- sort(unique(merged$true_effect[merged$true_effect > 0]))
  power <- data.table(effect_size = effect_sizes)
  for (nm in names(present_tests)) {
    cn <- present_tests[[nm]]
    power[[paste0("power_", nm)]] <- sapply(effect_sizes, function(eff) {
      d <- merged[true_effect == eff]
      if (nrow(d) == 0) return(NA_real_)
      mean(d[[cn]], na.rm = TRUE)
    })
  }

  structure(
    list(convergence = convergence, fpr = fpr, power = power),
    class = "reditR_validation"
  )
}


#' Print a reditR validation summary
#'
#' @param x An object of class \code{"reditR_validation"}.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#' @export
print.reditR_validation <- function(x, ...) {
  cat("reditR Validation Summary\n")
  cat("=========================\n")
  if (!is.na(x$convergence))
    cat(sprintf("  Convergence rate (GLMM):        %.1f%%\n", 100 * x$convergence))
  cat("  False-positive rate (null sites):\n")
  for (nm in names(x$fpr)) {
    cat(sprintf("    %-10s%.3f\n", paste0(nm, ":"),
               if (is.na(x$fpr[[nm]])) NaN else x$fpr[[nm]]))
  }
  cat("  Power by effect size:\n")
  print(x$power)
  invisible(x)
}
