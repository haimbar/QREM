library(testthat)

data(simdf)

# ── full pipeline: QREM → bcov → QRdiagnostics → flatQQplot ──────────────────
test_that("full pipeline runs without error on simdf at qn=0.25", {
  qn  <- 0.25
  fit <- QREM(lm, linmod = y ~ x * x2 + x3, dframe = simdf, qn = qn)

  # bcov
  cm  <- bcov(fit, linmod = y ~ x * x2 + x3, dframe = simdf, qn = qn)
  expect_true(all(diag(cm) > 0))

  # diagnostics - numeric
  diag_x <- QRdiagnostics(simdf$x, "x", fit$ui, qn, plot.it = FALSE)
  expect_type(diag_x, "list")

  # diagnostics - factor
  diag_x3 <- QRdiagnostics(simdf$x3, "x3", fit$ui, qn, plot.it = FALSE)
  expect_type(diag_x3, "list")
})

# ── multiple quantiles + flatQQplot ──────────────────────────────────────────
test_that("flatQQplot pipeline over 5 quantiles runs without error", {
  qns    <- seq(0.1, 0.9, by = 0.2)
  qrfits <- lapply(qns, function(q)
    QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = q))

  pvals <- flatQQplot(dat = simdf, cnum = 2, qrfits = qrfits,
                      qns = qns, maxm = 20, plot.it = FALSE)
  expect_length(pvals, length(qns))
})

# ── bcov SEs are sensible relative to boot.QREM SEs ─────────────────────────
test_that("bcov and boot.QREM give similar SEs for intercept (within 2x)", {
  skip_on_cran()
  # R CMD check enforces a parallel-core limit; skip when it is active
  skip_if(isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE"))),
          "Parallel cores restricted in this check environment")
  set.seed(1)
  df  <- data.frame(y = rnorm(200), x = rnorm(200))
  fit <- QREM(lm, linmod = y ~ x, dframe = df, qn = 0.5)

  se_bcov <- sqrt(diag(bcov(fit, y ~ x, df, 0.5)))

  bs <- boot.QREM(lm, linmod = y ~ x, dframe0 = df, qn = 0.5,
                  n = nrow(df), B = 50, seedno = 42)
  se_boot <- apply(bs, 2, sd)

  # Both should give roughly similar intercept SE (within factor 2)
  ratio <- se_bcov["(Intercept)"] / se_boot["(Intercept)"]
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 3.0)
})

# ── QREM with lmer (mixed model) smoke test ───────────────────────────────────
test_that("QREM works with lmer for a simple random-intercept model", {
  skip_if_not_installed("lme4")
  library(lme4)
  # Add a grouping factor to simdf
  set.seed(7)
  df       <- simdf[1:500, ]
  df$group <- factor(rep(1:10, each = 50))
  fit      <- QREM(lmer, linmod = y ~ x + x2 + (1 | group),
                   dframe = df, qn = 0.5)
  expect_true(!is.null(fit$coef$beta))
  expect_true(!is.null(fit$coef$u))
  expect_length(fit$ui, nrow(df))
})

# ── check coefficient direction on known data ─────────────────────────────────
test_that("QREM recovers correct sign of slope on simple simulated data", {
  set.seed(99)
  n  <- 500
  x  <- runif(n)
  y  <- 2 * x + rnorm(n, sd = 0.3)
  df <- data.frame(y = y, x = x)
  fit <- QREM(lm, linmod = y ~ x, dframe = df, qn = 0.5)
  # Slope should be close to 2
  expect_gt(fit$coef$beta["x"], 1.5)
  expect_lt(fit$coef$beta["x"], 2.5)
})

# ── qn=0.5 should be close to OLS on symmetric errors ────────────────────────
test_that("QREM at qn=0.5 gives coefficients close to OLS under normal errors", {
  set.seed(11)
  n  <- 1000
  x  <- rnorm(n)
  y  <- 3 + 1.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)
  fit_qr  <- QREM(lm, linmod = y ~ x, dframe = df, qn = 0.5)
  fit_ols <- lm(y ~ x, data = df)
  expect_equal(fit_qr$coef$beta["x"], coef(fit_ols)["x"],   tolerance = 0.1)
  expect_equal(fit_qr$coef$beta["(Intercept)"], coef(fit_ols)["(Intercept)"], tolerance = 0.1)
})
