library(testthat)

data(simdf)

# ── return type ───────────────────────────────────────────────────────────────
test_that("bcov returns a matrix", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  cm  <- bcov(fit, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_true(is.matrix(cm))
})

# ── dimensions match number of coefficients ───────────────────────────────────
test_that("bcov dimensions equal number of fixed-effect coefficients", {
  fit  <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  cm   <- bcov(fit, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  p    <- length(fit$coef$beta)
  expect_equal(dim(cm), c(p, p))
})

# ── symmetry ──────────────────────────────────────────────────────────────────
test_that("bcov returns a symmetric matrix", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  cm  <- bcov(fit, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_equal(cm, t(cm), tolerance = 1e-10)
})

# ── positive diagonal ─────────────────────────────────────────────────────────
test_that("bcov diagonal elements are positive", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  cm  <- bcov(fit, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  expect_true(all(diag(cm) > 0))
})

# ── SEs are finite and positive ───────────────────────────────────────────────
test_that("bcov standard errors are finite and positive", {
  fit <- QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  cm  <- bcov(fit, linmod = y ~ x + x2 + x3, dframe = simdf, qn = 0.5)
  se  <- sqrt(diag(cm))
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})

# ── varies sensibly with qn (extreme quantile → larger SE) ───────────────────
test_that("bcov SE is larger at qn=0.1 than at qn=0.5", {
  fit50 <- QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = 0.5)
  fit10 <- QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = 0.1)
  se50  <- sqrt(diag(bcov(fit50, y ~ x + x2, simdf, 0.5)))
  se10  <- sqrt(diag(bcov(fit10, y ~ x + x2, simdf, 0.1)))
  expect_gt(se10[1], se50[1])
})

# ── fix #5: largestNeg index bug — verify f(0) interpolation is correct ───────
test_that("bcov f(0) estimate is positive for well-behaved residuals", {
  # For a clean normal fit, f(0) of N(0, sigma) should be 1/(sigma*sqrt(2*pi))
  set.seed(42)
  df   <- data.frame(y = rnorm(500), x = rnorm(500))
  fit  <- QREM(lm, linmod = y ~ x, dframe = df, qn = 0.5)
  cm   <- bcov(fit, linmod = y ~ x, dframe = df, qn = 0.5)
  # SE of intercept should be finite and in a reasonable range
  se_int <- sqrt(diag(cm)["(Intercept)"])
  expect_true(is.finite(se_int))
  expect_gt(se_int, 0)
  expect_lt(se_int, 1)   # for n=500, SE << 1
})

# ── fix #6: near-zero dFp0 triggers a warning ─────────────────────────────────
test_that("bcov warns when estimated density at zero is near zero", {
  fit <- QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = 0.5)
  # Bimodal residuals with a gap at zero: bkde density at 0 will be tiny
  fit$ui <- c(rep(-5, length(fit$ui)/2), rep(5, length(fit$ui)/2))
  expect_warning(
    bcov(fit, linmod = y ~ x + x2, dframe = simdf, qn = 0.5),
    "near zero"
  )
})

# ── works across quantiles ────────────────────────────────────────────────────
test_that("bcov works for qn = 0.25, 0.5, 0.75", {
  for (qn in c(0.25, 0.5, 0.75)) {
    fit <- QREM(lm, linmod = y ~ x + x2, dframe = simdf, qn = qn)
    cm  <- bcov(fit, linmod = y ~ x + x2, dframe = simdf, qn = qn)
    expect_true(all(diag(cm) > 0),
                label = paste("positive diagonal at qn =", qn))
  }
})
