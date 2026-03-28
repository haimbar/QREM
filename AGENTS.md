# QREM – Agent Guide

This file helps AI agents (Claude Code, Copilot, Cursor, etc.) understand the
QREM package quickly so they can assist with development, testing, and
maintenance without needing to read every source file first.

---

## Package Overview

**QREM** fits quantile regression models via an EM algorithm.  The statistical
approach treats the problem as maximum likelihood under the Asymmetric Laplace
Distribution (ALD).  The EM iterations alternate between:

- **E-step** – compute per-observation weights as `1 / |residual|`
- **M-step** – fit a weighted least-squares model on a shifted response

The package extends ordinary quantile regression to support **random effects**
(`lmer`) and **variable selection** (`QREM_vs`).

Reference: <https://arxiv.org/abs/2104.08595>

---

## Repo Layout

```
R/QREM.R          All source code (single file, ~700 lines)
tests/testthat/   Test suite (79 tests across 5 files)
data/             simdf — built-in simulation dataset
doc/QREM.pdf      Full package manual
man/              Roxygen-generated .Rd files (do not edit manually)
DESCRIPTION       Package metadata, version, dependencies
NEWS.md           Changelog (prepend new sections at top)
```

---

## Exported Functions

### `QREM(func, linmod, dframe, qn, ...)`

Main fitting function.

| Argument | Description |
|----------|-------------|
| `func` | Fitting function for M-step: `lm`, `lmer`, or `gam` |
| `linmod` | Formula, e.g. `y ~ x + x2 * x3` |
| `dframe` | Data frame containing all formula variables |
| `qn` | Target quantile, in (0, 1) |
| `userwgts` | Optional column name/index for sampling weights |
| `err` | Initial log-likelihood diff to enter loop (default 10) |
| `maxit` | Max EM iterations (default 1000) |
| `tol` | Convergence threshold on log-likelihood change (default 0.001) |
| `maxInvLambda` | Cap on E-step weights (default 300) |

**Returns** a list: `coef`, `fitted.mod`, `empq`, `ui`, `weights`, `iter`, `err`

- `coef$beta` — fixed-effects coefficient vector (named)
- `coef$u` — random effects (lmer only)
- `fitted.mod` — the final M-step model object; pass to `summary()`, `aov()`, etc.
- `ui` — QR residuals (`y - fitted`); needed by `bcov` and `QRdiagnostics`
- `empq` — empirical fraction below the fit; should be ≈ `qn` at convergence

**Implementation note (lm fast path):** When `func = lm`, the inner loop solves
WLS via normal equations (`crossprod` + Cholesky) instead of calling `lm.wfit`.
This is ~1.5× faster per iteration.  A `tryCatch` falls back to `lm.wfit` on
Cholesky failure (near-singular design).

---

### `bcov(qremFit, linmod, dframe, qn, userwgts = NULL)`

Analytic covariance matrix for fixed-effects coefficients, using Bahadur's
representation:

```
Cov(β̂) = τ(1-τ) / f(0)² · (X'X)⁻¹
```

where `f(0)` is estimated by kernel density (Silverman bandwidth).

**Returns** a named `p × p` covariance matrix.  Take `sqrt(diag(result))` for
standard errors.  Named rows/columns match `coef(qremFit$fitted.mod)`.

**Prefer `bcov` over `boot.QREM` for fixed-effects models** (much faster).
For mixed models, use `boot.QREM`.

---

### `boot.QREM(func, linmod, dframe0, qn, n, B = 100, ...)`

Non-parametric bootstrap SE estimation.  Runs `B` QREM fits in parallel via
`parallel::parLapply`.

**Returns** a `B × p` matrix.  Get SEs with `apply(result, 2, sd)`.

Useful when `bcov` is unreliable (mixed models, very small samples).

---

### `QRdiagnostics(X, varname, u_i, qn, plot.it = TRUE, filename = NULL)`

Residual diagnostic for one predictor.

- **Numeric `X`**: produces a QQ-plot of binned residuals
- **Factor `X`**: computes the proportion below the regression line per level
  (spinogram)

**Returns** a list with `$dev` (a deviance-like scalar) and QQ coordinates
(`$x`, `$y`) for numeric predictors.

Returns `NULL` with a `warning()` for:
- constant residuals (zero variance after centering)
- unsupported predictor types (e.g., `character`)

**Typical usage:**
```r
fit <- QREM(lm, y ~ x + x2, dframe, qn = 0.5)
QRdiagnostics(dframe$x, "x", fit$ui, 0.5)
```

---

### `flatQQplot(dat, cnum = NULL, vname = NULL, qrfits, qns, maxm = 20, plot.it = TRUE)`

Heatmap diagnostic across multiple quantiles simultaneously.  Cell colour
represents standardised deviation of the observed proportion below the
regression line from the expected proportion `qn`, tested via `prop.test`.

- `dat` — data frame
- `cnum` — column index of the predictor to plot against
- `vname` — column name (alternative to `cnum`)
- `qrfits` — list of QREM fit objects (one per quantile)
- `qns` — vector of quantiles matching `qrfits`
- `maxm` — max number of bins/levels on the x-axis

**Returns** a list of `prop.test` result objects (one per quantile).

---

### `QREM_vs(inputData, ycol, Zcols, Xcols, qn, ...)`

Variable selection for high-dimensional quantile regression.  Uses the GAM
framework (alternating minimisation) with optional initialisation from the
`edgefinder` package.

**Returns** a list with selected variable indices and fitted model.

---

## Built-in Dataset

```r
data(simdf)
# 1000 rows × 5 columns: y, x, x2, x3 (factor), x4
# y is generated as a nonlinear function of x, x2, x3 with heteroscedastic noise
```

`simdf` is used throughout the test suite and in all `@examples`.

---

## Test Suite

```
tests/testthat/test-QREM.R           QREM return structure and convergence
tests/testthat/test-bcov.R           bcov covariance properties
tests/testthat/test-QRdiagnostics.R  QRdiagnostics numeric/factor/edge cases
tests/testthat/test-flatQQplot.R     flatQQplot output and edge cases
tests/testthat/test-integration.R    Full pipeline + coefficient recovery
```

Run tests from R:
```r
devtools::test()         # all tests
devtools::check()        # full R CMD check (includes tests, examples, docs)
```

Expected result: `[ FAIL 0 | WARN 1 | SKIP 0 | PASS 79 ]`
(The one warning is a harmless lme4 version mismatch on some systems.)

---

## Common Development Tasks

### Adding a new exported function

1. Add the function with Roxygen2 docs to `R/QREM.R`
2. Run `devtools::document()` to regenerate `NAMESPACE` and `man/`
3. Add tests in the appropriate `tests/testthat/test-*.R` file
4. Add an entry to `NEWS.md` under the current version section

### Bumping the version

1. Edit `Version:` in `DESCRIPTION`
2. Prepend a new `# QREM x.y.z` section to `NEWS.md`

### Checking for R CMD check issues

```r
devtools::check()
```

Common NOTEs to avoid:
- Use `inherits(x, "classname")` / `is.factor()` / `is.numeric()` instead of
  `class(x) == "classname"`
- All non-base package functions must be listed under `@importFrom` or in
  `DESCRIPTION`'s `Imports`

---

## Key Dependencies

| Package | Used for |
|---------|----------|
| `lme4` | Mixed-effects M-step (`lmer`), `getME`, `fixef`, `ranef` |
| `KernSmooth` | Kernel density estimate in `bcov` (`bkde`) |
| `gam` | GAM M-step in `QREM_vs` |
| `corpcor` | `make.positive.definite` fallback in `bcov` |
| `SEMMS` | Variable selection utilities in `QREM_vs` |
| `edgefinder` | Optional initialisation in `QREM_vs` |
| `Matrix` | Sparse diagonal matrix in `bcov` weighted path |
| `parallel` | Parallel bootstrap in `boot.QREM` |

---

## Gotchas

- **`chol2inv` does not preserve `dimnames`**.  After `chol2inv(chol(M))`, the
  result matrix has generic `[,1]`, `[,2]` column names.  If you need named
  subsetting (e.g. `diag(cm)["(Intercept)"]`), restore dimnames explicitly.
  See the `bcov` implementation for the pattern.

- **`tabulate` vs `table` for factor bins**.  `table` on a factor builds a
  full factor object internally; `tabulate(as.integer(x), nlevels(x))` is
  100× faster for large vectors.  Use `tabulate` wherever only bin counts are
  needed.

- **Hoist invariant work out of quantile loops**.  In `flatQQplot`, the `cut()`
  and bin-count steps do not depend on which quantile is being processed —
  compute them once before the loop.

- **The `lm` fast path uses normal equations**.  The general path (lmer, gam)
  still calls the fitting function directly in the loop and uses `logLik` for
  convergence.

- **`maxInvLambda`** caps E-step weights at 300.  If residuals are near zero,
  uncapped weights cause ill-conditioned WLS; don't remove this guard.
