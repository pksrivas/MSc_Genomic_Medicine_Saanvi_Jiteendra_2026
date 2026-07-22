#' Calculate per-site editing difference between conditions
#'
#' Computes the mean editing ratio per condition for each site and returns
#' the signed difference (case minus reference). Sites are ranked by absolute
#' editing difference. This is the R package equivalent of
#' \code{editing_difference.r} in the original pipeline.
#'
#' @param data_path Path to the filtered editing count file
#'   (\code{filtered_sites_clustered.txt} or similar). Required columns:
#'   \code{site}, \code{sample}, \code{edit_ratio}, \code{condition}.
#'   If \code{condition} is absent the metadata file is required via
#'   \code{meta_path}.
#' @param meta_path Path to the sample metadata file
#'   (\code{sample_metadata.txt}). Required columns: \code{sample},
#'   \code{condition}. Ignored if \code{condition} is already present in the
#'   data file.
#' @param out_path File path for the tab-delimited results table. Pass
#'   \code{NULL} (default) to skip writing.
#' @param reference_level Condition label for the reference group. Default
#'   \code{"control"}.
#' @param case_level Condition label for the case / disease group. Default
#'   \code{"diabetic"}.
#'
#' @return A \code{data.table} with columns \code{site},
#'   \code{<case_level>_mean}, \code{<reference_level>_mean},
#'   \code{editing_difference} (= case − reference), sorted by
#'   \code{abs(editing_difference)} descending.
#'
#' @examples
#' \dontrun{
#' eff <- editing_difference(
#'   data_path = "filtered_sites_clustered.txt",
#'   meta_path = "sample_metadata.txt"
#' )
#' head(eff)
#' }
#'
#' @export
editing_difference <- function(data_path,
                                meta_path      = NULL,
                                out_path       = NULL,
                                reference_level = "control",
                                case_level     = "diabetic") {

  # Load editing counts
  data <- fread(data_path)

  # Compute edit_ratio if absent
  if (!"edit_ratio" %in% names(data))
    data[, edit_ratio := edited / total]

  # Merge with metadata if condition column is absent
  if (!"condition" %in% names(data)) {
    if (is.null(meta_path))
      stop("'condition' column not found in data file and no meta_path supplied.")
    meta <- fread(meta_path)
    data <- merge(data, meta, by = "sample")
  }

  # Column names derived from the supplied condition labels
  case_col <- paste0(case_level, "_mean")
  ref_col  <- paste0(reference_level, "_mean")

  # Calculate mean editing per condition
  editing_summary <- data[, .(
    case_mean = mean(edit_ratio[condition == case_level],      na.rm = TRUE),
    ref_mean  = mean(edit_ratio[condition == reference_level], na.rm = TRUE)
  ), by = site]

  setnames(editing_summary, c("case_mean", "ref_mean"), c(case_col, ref_col))

  # Calculate difference
  editing_summary$editing_difference <-
    editing_summary[[case_col]] - editing_summary[[ref_col]]

  # Rank sites by absolute editing change
  editing_summary <- editing_summary[order(-abs(editing_difference))]

  if (!is.null(out_path))
    fwrite(editing_summary, out_path, sep = "\t")

  cat("Finished calculating editing differences\n")
  editing_summary
}
