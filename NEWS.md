# QREM 0.2.0

## Bug fixes

* `QREM()`: `modelCoefs` and `ui` are now initialised from the pre-loop fit,
  preventing a crash when `err = 0` is passed and the EM loop is skipped.
* `QREM()`: Added an explicit `stop()` if the response variable named in
  `linmod` is not found in `dframe`, replacing a silent wrong-data error.
* `boot.QREM()`: `sampleFrom` is now included in `clusterExport()`, fixing a
  crash that occurred whenever `sampleFrom` was non-`NULL` and parallel
  workers tried to access it.
* `boot.QREM()`: `detectCores() - 1` is now guarded with `na.rm = TRUE`,
  preventing a crash on systems (including some Windows configurations) where
  `detectCores()` returns `NA`.
* `boot.QREM()`: Bootstrap replicates whose coefficient vector has a different
  length than the initial fit (due to rank deficiency in a resample) are now
  silently dropped before `rbind()`, preventing a non-conformable-dimensions
  error.
* `bcov()`: Fixed an indexing bug in the kernel-density interpolation of
  `f(0)`: `which.max()` was applied to a filtered sub-vector, producing an
  index into that sub-vector that was then used as an index into the full KDE
  grid, yielding the wrong interpolation points. The fix uses `which()` to
  obtain the correct position in the full grid.
* `bcov()`: Added a `warning()` when the estimated density at zero (`dFp0`) is
  near machine epsilon, alerting users to unreliable covariance estimates
  instead of silently returning enormous values.
* `bcov()`: When a mixed model (`lmerMod`) has only one random-effects column,
  `ZZ[, -ncol(ZZ)]` previously dropped all columns. The design matrix now
  falls back to `XX` alone in that edge case.
* `QRdiagnostics()`: Added a guard for zero-variance residuals (all identical
  after centering); the function now returns `NULL` with a warning instead of
  crashing inside `bkde()` with a zero-bandwidth error.
* `QRdiagnostics()`: Added an early return with `warning()` and `dev.off()` for
  unsupported predictor types (e.g., `character`), closing any open PDF device
  and preventing a silent `NULL` return without explanation.
* `flatQQplot()`: The factor-column branch now applies the same one-row table
  guard that the numeric branch already had, preventing a subscript
  out-of-bounds crash when all observations fall on one side of the regression
  line for a given quantile.
* `flatQQplot()`: The sample-size guard now catches `m == 0` in addition to
  `m == 1`, preventing a `'probs' outside [0,1]` crash when the dataset is
  very small relative to the number of quantiles requested. Also fixed a typo
  (`warnings()` → `warning()`).
* `flatQQplot()`: `NA` values returned by `cut()` (when a standardised cell
  count falls outside the colour-scale range) are now mapped to the most
  extreme colour rather than causing an `NA`-subscript error in `rect()`.
* `QREM_vs()`: `fittedVSnew` is now initialised to `NULL` before the main
  loop, preventing a crash when the convergence criterion is met on the very
  first iteration.
* `QREM_vs()`: The initialisation loop over groups of 5 predictors now uses
  `ceiling(K / m)` with `min(i * m, K)` as the upper index, ensuring all `K`
  predictors are scored when `K` is not divisible by 5.
* `QREM_vs()`: Standard errors in the initial z-score calculation are now
  computed as `sqrt(pmax(diag(bcov(...)), 0))`, preventing `NaN` propagation
  from any numerically negative diagonal entry.

## Improvements

* All `class(x) == "..."` comparisons replaced with `inherits()`,
  `is.factor()`, `is.numeric()`, and `is.data.frame()` throughout, resolving
  `R CMD check` NOTEs and correctly treating integer vectors as numeric in
  `QRdiagnostics()`.
* `QREM()`, `bcov()`: Documentation substantially expanded — the ALD/EM
  formulation and Bahadur's representation formula are now stated explicitly,
  parameter descriptions clarified (e.g., `userwgts` now states column name
  or index; `maxInvLambda` explains the numerical-stability rationale), and
  `@seealso` cross-references added across functions.

## Testing

* Added a `tests/testthat/` suite (79 tests across 5 files) covering return
  structure, convergence, coefficient recovery, covariance matrix properties,
  diagnostic plots, heatmap output, and the full QREM → `bcov` →
  `QRdiagnostics` → `flatQQplot` pipeline.

## Package infrastructure

* Added `Suggests: testthat (>= 3.0.0)` to `DESCRIPTION`.
* Updated `RoxygenNote` to 7.3.2.
* Added `^doc$`, `^README\.md$`, `^\.DS_Store$`, and `^\.gitignore$` to
  `.Rbuildignore`, removing a spurious "non-standard top-level directory" NOTE
  from `R CMD check`.

---

# QREM 0.1.9

* Added `flatQQplot()`: a heatmap-based diagnostic showing goodness of fit
  across multiple quantiles simultaneously, with cell colours representing
  standardised deviations from the expected proportion and significance tested
  via `prop.test()`.
* Updated documentation.

# QREM 0.1.8

* Fixed initialisation of `QREM_vs()` when `initWithEdgeFinder = TRUE`.

# QREM 0.1.7

* Fixed a bug in `QREM_vs()` in the `edgefinder`-initialised path.

# QREM 0.1.0

* First public release on GitHub.
