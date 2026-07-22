#' Filter editing sites by coverage, read count, and cluster proximity
#'
#' Reads a merged editing count table, applies per-row QC filters, and then
#' identifies genomically clustered sites — sites within \code{cluster_window}
#' bp of a neighbouring site on the same chromosome. This is the R package
#' equivalent of \code{filtered_sites.R} in the original pipeline.
#'
#' The HPC/SPRINT helper script that produces the input file is shipped at
#' \code{system.file("scripts", "extracting_read_counts.sh", package = "reditR")}.
#'
#' @param data_path Path to the merged per-sample editing count file
#'   (\code{all_samples_editing.txt}). Required columns: \code{site},
#'   \code{sample}, \code{edited}, \code{total}.
#' @param out_dir Directory where \code{filtered_sites_all.txt} and
#'   \code{filtered_sites_clustered.txt} will be written. Pass \code{NULL}
#'   (default) to skip writing.
#' @param min_coverage Minimum total read coverage per observation. Default 10.
#' @param min_edited Minimum edited read count per observation. Default 2.
#' @param min_groups Minimum number of distinct samples a site must appear in
#'   after per-row filtering. Default 2.
#' @param cluster_window Maximum distance in bp between two sites on the same
#'   chromosome for them to be considered clustered. Default 50.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{all}}{All QC-passing sites.}
#'     \item{\code{clustered}}{Subset of sites within a genomic cluster.}
#'   }
#'
#' @examples
#' \dontrun{
#' res   <- filter_editing_sites("all_samples_editing.txt", out_dir = "results/")
#' sites <- res$clustered
#' }
#'
#' @export
filter_editing_sites <- function(data_path,
                                  out_dir        = NULL,
                                  min_coverage   = 10,
                                  min_edited     = 2,
                                  min_groups     = 2,
                                  cluster_window = 50) {

  # Load editing data
  data <- read.table(data_path, header = TRUE)

  # Calculate editing ratio
  data <- data %>%
    mutate(edit_ratio = edited / total)

  # Filter: coverage, edited reads, site in >= min_groups samples
  data_filtered <- data %>%
    filter(total >= min_coverage, edited >= min_edited) %>%
    group_by(site) %>%
    filter(n_distinct(sample) >= min_groups) %>%
    ungroup()

  # Prepare coordinates
  data_filtered <- data_filtered %>%
    separate(site, into = c("chr", "pos"), sep = ":")

  data_filtered$pos <- as.numeric(data_filtered$pos)

  # Identify clustered genomic sites
  unique_sites <- data_filtered %>%
    distinct(chr, pos)

  clustered_positions <- unique_sites %>%
    group_by(chr) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(
      dist_prev = pos - lag(pos),
      dist_next = lead(pos) - pos
    ) %>%
    filter(dist_prev <= cluster_window | dist_next <= cluster_window) %>%
    ungroup()

  # Keep rows belonging to clustered sites
  clustered_sites <- data_filtered %>%
    semi_join(clustered_positions, by = c("chr", "pos"))

  # Reconstruct site column
  data_filtered <- data_filtered %>%
    mutate(site = paste(chr, pos, sep = ":")) %>%
    select(site, sample, edited, total, edit_ratio, chr, pos)

  clustered_sites <- clustered_sites %>%
    mutate(site = paste(chr, pos, sep = ":")) %>%
    select(site, sample, edited, total, edit_ratio, chr, pos)

  # Save results
  if (!is.null(out_dir)) {
    write.table(data_filtered,
                file.path(out_dir, "filtered_sites_all.txt"),
                sep       = "\t",
                row.names = FALSE,
                quote     = FALSE)

    write.table(clustered_sites,
                file.path(out_dir, "filtered_sites_clustered.txt"),
                sep       = "\t",
                row.names = FALSE,
                quote     = FALSE)
  }

  invisible(list(all = data_filtered, clustered = clustered_sites))
}
