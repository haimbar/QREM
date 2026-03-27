library(testthat)

data(simdf)

fit_med <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)

# ── numeric predictor ─────────────────────────────────────────────────────────
test_that("QRdiagnostics returns qqp list with $dev for numeric predictor", {
  result <- QRdiagnostics(simdf$x, "x", fit_med$ui, 0.5, plot.it = FALSE)
  expect_type(result, "list")
  expect_true(!is.null(result$dev))
  expect_true(is.finite(result$dev))
})

test_that("QRdiagnostics qqp contains x and y quantile vectors", {
  result <- QRdiagnostics(simdf$x, "x", fit_med$ui, 0.5, plot.it = FALSE)
  expect_true(!is.null(result$x))
  expect_true(!is.null(result$y))
})

# ── factor predictor ──────────────────────────────────────────────────────────
test_that("QRdiagnostics returns list with $dev for factor predictor", {
  result <- QRdiagnostics(simdf$x3, "x3", fit_med$ui, 0.5, plot.it = FALSE)
  expect_type(result, "list")
  expect_true(!is.null(result$dev))
  expect_true(is.finite(result$dev))
})

test_that("QRdiagnostics factor result has one entry per level plus $dev", {
  result <- QRdiagnostics(simdf$x3, "x3", fit_med$ui, 0.5, plot.it = FALSE)
  n_levels <- nlevels(simdf$x3)
  expect_equal(length(result) - 1L, n_levels)  # $dev is the extra element
})

# ── fix #8: constant residuals (sd=0 after centering) warns and returns NULL ──
test_that("QRdiagnostics warns and returns NULL for constant residuals", {
  ui_const <- rep(1.5, nrow(simdf))
  expect_warning(
    result <- QRdiagnostics(simdf$x, "x", ui_const, 0.5, plot.it = FALSE),
    "identical"
  )
  expect_null(result)
})

# ── fix #4: unsupported type warns and returns NULL ───────────────────────────
test_that("QRdiagnostics treats integer X as numeric (QQ-plot path)", {
  # is.numeric(integer) is TRUE in R, so integers correctly go through the QQ-plot path
  X_int <- as.integer(seq_len(nrow(simdf)))
  result <- QRdiagnostics(X_int, "int_var", fit_med$ui, 0.5, plot.it = FALSE)
  expect_type(result, "list")
  expect_true(!is.null(result$dev))
})

test_that("QRdiagnostics warns and returns NULL for character X", {
  X_chr <- as.character(simdf$x3)
  expect_warning(
    result <- QRdiagnostics(X_chr, "chr_var", fit_med$ui, 0.5, plot.it = FALSE),
    "numeric or factor"
  )
  expect_null(result)
})

# ── deviance is non-negative ──────────────────────────────────────────────────
test_that("QRdiagnostics deviance is non-negative", {
  result <- QRdiagnostics(simdf$x, "x", fit_med$ui, 0.5, plot.it = FALSE)
  expect_gte(result$dev, 0)
})

# ── works across quantiles ────────────────────────────────────────────────────
test_that("QRdiagnostics works for qn=0.1 and qn=0.9", {
  for (qn in c(0.1, 0.9)) {
    fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = qn)
    result <- QRdiagnostics(simdf$x, "x", fit$ui, qn, plot.it = FALSE)
    expect_type(result, "list")
  }
})
