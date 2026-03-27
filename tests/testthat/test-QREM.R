library(testthat)

# ── helpers ──────────────────────────────────────────────────────────────────
data(simdf)

# ── return structure ──────────────────────────────────────────────────────────
test_that("QREM returns a list with all expected fields", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_type(fit, "list")
  expect_named(fit, c("coef", "fitted.mod", "empq", "ui", "weights", "iter", "err"),
               ignore.order = TRUE)
})

test_that("QREM coef list has beta for lm fit", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_true(!is.null(fit$coef$beta))
  expect_named(fit$coef, "beta")
})

# ── empirical quantile ────────────────────────────────────────────────────────
test_that("QREM empq is close to qn at convergence", {
  for (qn in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = qn)
    expect_equal(fit$empq, qn, tolerance = 0.05,
                 label = paste("empq for qn =", qn))
  }
})

# ── convergence ───────────────────────────────────────────────────────────────
test_that("QREM converges within maxit iterations", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_lt(fit$iter, 1000)
  expect_lt(fit$err, 0.001)
})

# ── residuals length ──────────────────────────────────────────────────────────
test_that("QREM residuals ui has same length as data", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_length(fit$ui, nrow(simdf))
})

# ── weights are positive and bounded ─────────────────────────────────────────
test_that("QREM weights are positive and at most maxInvLambda", {
  fit <- QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = 0.5,
              maxInvLambda = 300)
  expect_true(all(fit$weights > 0))
  expect_true(all(fit$weights <= 300))
})

# ── interaction terms ─────────────────────────────────────────────────────────
test_that("QREM works with interaction formula", {
  fit <- QREM(lm, linmod = y ~ x * x2 + x3, dframe = simdf, qn = 0.2)
  expect_s3_class(fit$fitted.mod, "lm")
})

# ── fix #16: missing response column stops with a clear message ───────────────
test_that("QREM stops if response column not in dframe", {
  expect_error(
    QREM(lm, linmod = z ~ x + x2, dframe = simdf, qn = 0.5),
    "not found in dframe"
  )
})

# ── fix #10: err=0 does not crash (modelCoefs/ui initialized pre-loop) ────────
test_that("QREM does not crash when err=0 (loop skipped)", {
  # With err=0 the while loop condition (err > tol) is FALSE immediately,
  # so the pre-loop initialisation must supply modelCoefs and ui.
  expect_no_error(
    QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = 0.5, err = 0)
  )
})

# ── fitted.mod is a usable lm object ─────────────────────────────────────────
test_that("QREM fitted.mod supports summary and aov", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_no_error(summary(fit$fitted.mod))
  expect_no_error(summary(aov(fit$fitted.mod)))
})

# ── small dataset smoke test ──────────────────────────────────────────────────
test_that("QREM works on small simulated dataset", {
  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  fit <- QREM(lm, linmod = y ~ x, dframe = df, qn = 0.5)
  expect_type(fit, "list")
  expect_length(fit$ui, 50)
})
