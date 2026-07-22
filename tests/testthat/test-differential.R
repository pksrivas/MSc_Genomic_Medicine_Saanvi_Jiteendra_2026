library(testthat)
library(reditR)
library(data.table)

editing_file  <- system.file("extdata", "example_editing.txt",  package = "reditR")
metadata_file <- system.file("extdata", "example_metadata.txt", package = "reditR")

test_that("differential_editing runs all three tests by default with expected column schema", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  expect_s3_class(res, "data.table")
  expected <- c("site",
               "glmm_pvalue",   "GLMM_FDR",   "GLMM_sig",
               "fisher_pvalue", "Fisher_FDR", "Fisher_sig",
               "wilcox_pvalue", "Wilcox_FDR", "Wilcox_sig")
  expect_true(all(expected %in% names(res)))
  expect_false("DRE" %in% names(res))
})

test_that("test= restricts which columns are computed", {
  res_glmm <- differential_editing(editing_file, metadata_file,
                                   test = "glmm", verbose = FALSE)
  expect_true(all(c("glmm_pvalue", "GLMM_sig") %in% names(res_glmm)))
  expect_false(any(c("fisher_pvalue", "wilcox_pvalue") %in% names(res_glmm)))

  res_wf <- differential_editing(editing_file, metadata_file,
                                 test = c("wilcoxon", "fisher"), verbose = FALSE)
  expect_true(all(c("wilcox_pvalue", "fisher_pvalue") %in% names(res_wf)))
  expect_false("glmm_pvalue" %in% names(res_wf))
})

test_that("chr1_100 (clear effect) is GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)
  expect_true(res[site == "chr1_100", GLMM_sig])
})

test_that("chr1_200 (no effect) is not GLMM significant", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)
  expect_false(res[site == "chr1_200", GLMM_sig])
})

test_that("each requested test is independently FDR-corrected (no gating, no combined verdict)", {
  res <- differential_editing(editing_file, metadata_file, verbose = FALSE)

  # Wilcoxon and Fisher must be computed for every site with >=2 conditions,
  # not just GLMM-significant ones (i.e. not gated on GLMM_sig).
  expect_equal(sum(!is.na(res$wilcox_pvalue)), nrow(res))
  expect_equal(sum(!is.na(res$fisher_pvalue)), nrow(res))

  # FDR correction is independent per test.
  expect_equal(res$GLMM_FDR,   p.adjust(res$glmm_pvalue,   method = "BH"))
  expect_equal(res$Fisher_FDR, p.adjust(res$fisher_pvalue, method = "BH"))
  expect_equal(res$Wilcox_FDR, p.adjust(res$wilcox_pvalue, method = "BH"))
})

test_that("results are written to out_path when provided", {
  out <- tempfile(fileext = ".txt")
  on.exit(unlink(out))

  differential_editing(editing_file, metadata_file, out_path = out, verbose = FALSE)

  expect_true(file.exists(out))
  written <- fread(out)
  expect_equal(nrow(written), 2L)
})

test_that("summary file is written with expected rows when summary_path is provided", {
  sumf <- tempfile(fileext = ".txt")
  on.exit(unlink(sumf))

  differential_editing(editing_file, metadata_file, summary_path = sumf, verbose = FALSE)

  expect_true(file.exists(sumf))
  sumdt <- fread(sumf)
  expect_true(all(c("metric", "value") %in% names(sumdt)))
  expect_true("GLMM_sig" %in% sumdt$metric)
  expect_true("Fisher_sig" %in% sumdt$metric)
  expect_true("Wilcox_sig" %in% sumdt$metric)
  expect_false("monte_carlo_pvalue" %in% sumdt$metric)
})
