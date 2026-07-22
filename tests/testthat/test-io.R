library(testthat)
library(reditR)
library(data.table)

editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")

test_that("read_editing_table reads correctly and computes edit_ratio if absent", {
  dt <- read_editing_table(editing_file)

  expect_s3_class(dt, "data.table")
  expect_true(all(c("site", "sample", "edited", "total", "edit_ratio") %in% names(dt)))
  expect_equal(dt[site == "chr1_100" & sample == "ctrl_1", edit_ratio], 2 / 30)
})

test_that("read_editing_table errors on missing required columns", {
  f <- tempfile(fileext = ".txt")
  fwrite(data.table(site = "s1", edited = 1L), f, sep = "\t")
  on.exit(unlink(f))

  expect_error(read_editing_table(f), "Missing required columns")
})

test_that("read_editing_table errors on nonexistent file", {
  expect_error(read_editing_table("/nonexistent/path.txt"), "File not found")
})

test_that("read_metadata reads correctly", {
  mt <- read_metadata(metadata_file)

  expect_s3_class(mt, "data.table")
  expect_true(all(c("sample", "condition") %in% names(mt)))
  expect_equal(nrow(mt), 6L)
})

test_that("read_metadata with reference_level sets factor levels correctly", {
  mt <- read_metadata(metadata_file, reference_level = "control")

  expect_true(is.factor(mt$condition))
  expect_equal(levels(mt$condition)[1], "control")
})

test_that("read_metadata errors when reference_level not in conditions", {
  expect_error(read_metadata(metadata_file, reference_level = "missing"), "not found")
})
