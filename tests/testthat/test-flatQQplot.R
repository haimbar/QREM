library(testthat)

data(simdf)

qns <- c(0.25, 0.5, 0.75)
qrfits <- lapply(qns, function(q)
  QREM(lm, linmod = y ~ x + x2 + x3, dframe = simdf, qn = q))

# ── return type ───────────────────────────────────────────────────────────────
test_that("flatQQplot returns a list of prop.test results for numeric column", {
  result <- flatQQplot(dat = simdf, cnum = 2, qrfits = qrfits, qns = qns,
                       maxm = 10, plot.it = FALSE)
  expect_type(result, "list")
  expect_length(result, length(qns))
})

test_that("flatQQplot returns a list of prop.test results for factor column", {
  result <- flatQQplot(dat = simdf, cnum = 4, qrfits = qrfits, qns = qns,
                       maxm = 10, plot.it = FALSE)
  expect_type(result, "list")
  expect_length(result, length(qns))
})

# ── vname works as alternative to cnum ───────────────────────────────────────
test_that("flatQQplot accepts vname instead of cnum", {
  result <- flatQQplot(dat = simdf, vname = "x", qrfits = qrfits, qns = qns,
                       maxm = 10, plot.it = FALSE)
  expect_type(result, "list")
})

# ── missing both cnum and vname returns FALSE with warning ────────────────────
test_that("flatQQplot warns and returns FALSE when cnum and vname both NULL", {
  expect_warning(
    result <- flatQQplot(dat = simdf, qrfits = qrfits, qns = qns,
                         plot.it = FALSE),
    "Must specify"
  )
  expect_false(result)
})

# ── fix #3: factor column where all obs fall on one side doesn't crash ────────
test_that("flatQQplot does not crash for factor column with one-sided residuals", {
  fits_biased <- qrfits
  fits_biased[[1]]$ui <- -abs(fits_biased[[1]]$ui) - 1
  expect_no_error(
    flatQQplot(dat = simdf, cnum = 4, qrfits = fits_biased, qns = qns,
               maxm = 10, plot.it = FALSE)
  )
})

# ── p-values are in [0, 1] ────────────────────────────────────────────────────
test_that("flatQQplot prop.test p-values are in [0,1]", {
  result <- flatQQplot(dat = simdf, cnum = 2, qrfits = qrfits, qns = qns,
                       maxm = 10, plot.it = FALSE)
  pvals <- sapply(result, function(r) r$p.value)
  expect_true(all(pvals >= 0 & pvals <= 1))
})

# ── sample-size guard: m==1 returns FALSE ────────────────────────────────────
test_that("flatQQplot returns FALSE when sample size too small for k quantiles", {
  # Use only numeric predictors on a tiny dataset to avoid dropped factor levels
  set.seed(1)
  tiny    <- data.frame(y = rnorm(15), x = rnorm(15))
  qns9    <- seq(0.1, 0.9, by = 0.1)
  qrfits9 <- lapply(qns9, function(q)
    QREM(lm, y ~ x, dframe = tiny, qn = q))
  expect_false(suppressWarnings(
    flatQQplot(dat = tiny, cnum = 2, qrfits = qrfits9, qns = qns9,
               maxm = 30, plot.it = FALSE)
  ))
})

# ── each element of pvals is a htest object ──────────────────────────────────
test_that("flatQQplot pvals list elements are htest objects", {
  result <- flatQQplot(dat = simdf, cnum = 2, qrfits = qrfits, qns = qns,
                       maxm = 10, plot.it = FALSE)
  for (r in result) {
    expect_s3_class(r, "htest")
  }
})
