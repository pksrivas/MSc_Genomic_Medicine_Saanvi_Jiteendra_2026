library(testthat)
library(reditR)
library(data.table)

editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")

test_that("editing_difference returns expected schema", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_s3_class(eff, "data.table")
  expect_true(all(c("site", "diabetic_mean", "control_mean",
                     "editing_difference") %in% names(eff)))
})

test_that("output column names reflect the supplied condition labels", {
  eff <- editing_difference(editing_file, meta_path = metadata_file,
                             case_level = "diabetic", reference_level = "control")
  expect_true("diabetic_mean" %in% names(eff))
  expect_true("control_mean"  %in% names(eff))
})

test_that("editing difference sign is positive for chr1_100 (diabetic > control)", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_gt(eff[site == "chr1_100", editing_difference], 0)
})

test_that("chr1_200 editing difference is near zero", {
  eff <- editing_difference(editing_file, meta_path = metadata_file)

  expect_lt(abs(eff[site == "chr1_200", editing_difference]), 0.05)
})

test_that("sites are sorted by decreasing absolute editing difference", {
  eff   <- editing_difference(editing_file, meta_path = metadata_file)
  diffs <- abs(eff$editing_difference)

  expect_true(all(diff(diffs) <= 0))
})

test_that("results are written to out_path when provided", {
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out))

  editing_difference(editing_file, meta_path = metadata_file, out_path = out)

  expect_true(file.exists(out))
})
