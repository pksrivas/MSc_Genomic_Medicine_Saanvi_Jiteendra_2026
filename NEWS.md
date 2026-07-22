# reditR 0.2.0

* **Breaking change:** `differential_editing()` redesigned around three
  independent, user-selectable tests rather than a single GLMM-primary /
  Wilcoxon-robustness-flag framework:
  - New `test` argument (any subset of `c("glmm", "fisher", "wilcoxon")`,
    default all three). Only the requested test(s) are computed.
  - All requested tests now run on **every** site, independently of one
    another -- Wilcoxon is no longer gated on GLMM significance.
  - Each requested test gets its own independent Benjamini-Hochberg FDR
    correction (`GLMM_FDR`/`GLMM_sig`, `Fisher_FDR`/`Fisher_sig`,
    `Wilcox_FDR`/`Wilcox_sig`).
  - Fisher's exact test is now a first-class option (previously absent).
  - Removed: the global Monte Carlo permutation overlap test (and its
    `n_perm`/`seed` arguments and `monte_carlo_pvalue` attribute), the
    combined `DRE` verdict column, and the raw-p-value Wilcoxon
    "robustness flag" convention. Requesting more than one test lets you
    compare columns directly, or build your own combined criterion.
  - Guidance on which test to pick for a given experimental design (based
    on null-calibration simulations across three independent datasets) is
    now documented in `?differential_editing` and the package vignette.
* `validate_against_truth()`/`print.reditR_validation()` updated to report
  false-positive rate and power independently per test present in the
  `differential_editing()` output, instead of hard-coding GLMM + DRE.
* Removed the `permutation_pvalue()` export: it was declared in `NAMESPACE`
  and documented, but never actually implemented in this package -- calling
  it would have errored. No functional change (it was already unusable).

# reditR 0.1.0

* Initial release.
* Core functions: `read_editing_table()`, `read_metadata()`,
  `filter_editing_sites()`, `differential_editing()`,
  `editing_difference()`, `simulate_editing_data()`,
  `validate_against_truth()`.
* Implements the GLMM + per-site robustness-flag framework for
  small-sample differential editing analysis.
