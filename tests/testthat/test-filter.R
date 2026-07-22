library(testthat)
library(reditR)
library(data.table)

# Write a temp editing file and run filter_editing_sites() on it
make_temp_editing_file <- function(rows) {
  f <- tempfile(fileext = ".txt")
  write.table(rows, f, sep = "\t", row.names = FALSE, quote = FALSE)
  f
}

base_rows <- data.frame(
  site   = c("chr1:100", "chr1:100", "chr1:100",
             "chr1:150", "chr1:150", "chr1:150",
             "chr2:900"),
  sample = paste0("S", 1:7),
  edited = c(5L,  3L, 4L,  2L, 1L, 3L, 1L),
  total  = c(30L, 25L, 20L, 15L, 12L, 18L, 5L)
)

test_that("coverage filter removes rows below min_coverage", {
  f    <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res  <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10, min_edited = 1, min_groups = 1)

  expect_true(all(res$all$total >= 10))
})

test_that("min_edited filter removes rows with too few edited reads", {
  f    <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res  <- filter_editing_sites(f, out_dir = NULL, min_coverage = 1, min_edited = 3, min_groups = 1)

  expect_true(all(res$all$edited >= 3))
})

test_that("min_groups filter removes underpresent sites", {
  single <- data.frame(site = "chr1:100", sample = "S1", edited = 5L, total = 30L)
  f      <- make_temp_editing_file(single)
  on.exit(unlink(f))
  res    <- filter_editing_sites(f, out_dir = NULL, min_coverage = 1, min_edited = 1, min_groups = 2)

  expect_equal(nrow(res$all), 0L)
})

test_that("clustered sites are a subset of all filtered sites", {
  f   <- make_temp_editing_file(base_rows)
  on.exit(unlink(f))
  res <- filter_editing_sites(f, out_dir = NULL, min_coverage = 10,
                               min_edited = 1, min_groups = 1, cluster_window = 100)

  expect_true(all(res$clustered$site %in% res$all$site))
})

test_that("output is written when out_dir is provided", {
  f      <- make_temp_editing_file(base_rows)
  tmpdir <- tempdir()
  on.exit(unlink(f))

  filter_editing_sites(f, out_dir = tmpdir, min_coverage = 10, min_edited = 1, min_groups = 1)

  expect_true(file.exists(file.path(tmpdir, "filtered_sites_all.txt")))
  expect_true(file.exists(file.path(tmpdir, "filtered_sites_clustered.txt")))
})
