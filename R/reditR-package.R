#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import tidyr
#' @importFrom data.table ":=" ".N" as.data.table data.table fread fwrite is.data.table rbindlist setnames uniqueN
#' @importFrom lme4 glmer
#' @importFrom stats binomial fisher.test p.adjust plogis qlogis rbinom rnbinom rnorm wilcox.test
#' @importFrom utils read.table write.table
## usethis namespace: end
NULL

utils::globalVariables(c(
  ".",
  "site", "sample", "edited", "total", "edit_ratio", "unedited", "condition",
  "chr", "pos", "dist_prev", "dist_next",
  "GLMM_FDR", "GLMM_sig", "Fisher_FDR", "Fisher_sig", "Wilcox_FDR", "Wilcox_sig",
  "glmm_pvalue", "fisher_pvalue", "wilcox_pvalue",
  "diabetic_mean", "control_mean", "editing_difference",
  "reference_mean", "other_mean",
  "true_effect", "N"
))
