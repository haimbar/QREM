#' Fitting a quantile regression model via the EM algorithm.
#'
#' Fits a quantile regression model by treating it as a maximum likelihood
#' problem under the Asymmetric Laplace Distribution (ALD). The EM algorithm
#' iterates between an E-step, which computes observation weights as
#' \eqn{1/|u_i|} (inversely proportional to the absolute residuals), and an
#' M-step, which fits a weighted least squares model on a shifted response.
#' Convergence is assessed by the change in log-likelihood between iterations.
#' @param func The fitting function to use in the M-step. Typically \code{lm},
#'   \code{\link[lme4]{lmer}}, or \code{\link[gam]{gam}}. Any function that
#'   accepts \code{formula}, \code{data}, and \code{weights} arguments and
#'   returns an object with a \code{\link{logLik}} method is supported.
#' @param linmod A formula specifying the linear model for the M-step.
#' @param dframe A data frame containing all variables referenced in
#'   \code{linmod}.
#' @param qn The target quantile. Must be in (0, 1).
#' @param userwgts A column name or index in \code{dframe} for user-provided
#'   sampling weights (optional. Default=NULL). When supplied, EM weights are
#'   multiplied by these sampling weights at each iteration.
#' @param ... Additional arguments passed to \code{func} (except
#'   \code{formula}, \code{data}, and \code{weights}, which are set
#'   internally). For \code{gam}, the \code{family} argument must be supplied
#'   here.
#' @param err Initial log-likelihood difference used to enter the EM loop
#'   (default=10). Should be set larger than \code{tol} to ensure at least one
#'   iteration runs.
#' @param maxit Maximum number of EM iterations (default=1000).
#' @param tol Convergence threshold on the absolute change in log-likelihood
#'   between iterations (default=0.001).
#' @param maxInvLambda Upper bound on the E-step weights \eqn{1/|u_i|}
#'   (default=300). Prevents numerical instability when residuals are close to
#'   zero. Increase cautiously; very large values can cause ill-conditioned WLS
#'   fits.
#' @return A list containing the following
#' \itemize{
#'   \item \code{coef} A list of estimated coefficients. For \code{lm} and
#'     \code{gam}: \code{list(beta = coef vector)}. For \code{lmerMod}:
#'     \code{list(beta = fixed effects, u = random effects)}.
#'   \item \code{fitted.mod} The model object returned by \code{func} in the
#'     final EM iteration. Can be passed to \code{summary()}, \code{aov()},
#'     etc.
#'   \item \code{empq} The proportion of observations below the fitted
#'     regression line. Should be close to \code{qn} at convergence.
#'   \item \code{ui} The quantile regression residuals
#'     (\eqn{y_i - \hat{y}_i}) from the final iteration.
#'   \item \code{weights} The final E-step weights (\eqn{1/|u_i|}, capped at
#'     \code{maxInvLambda}). Passed to \code{\link{bcov}} for covariance
#'     estimation.
#'   \item \code{iter} The number of EM iterations performed.
#'   \item \code{err} The final log-likelihood change (convergence criterion).
#'     Values above \code{tol} indicate the algorithm did not converge within
#'     \code{maxit} iterations.
#' }
#' @seealso \code{\link{bcov}}, \code{\link{boot.QREM}},
#'   \code{\link{QRdiagnostics}}
#' @importFrom lme4 lmer fixef ranef
#' @export
#' @examples
#' data(simdf)
#' qremFit <-  QREM(lm, linmod=y~x*x2 +x3, dframe=simdf, qn=0.2)
#' summary(aov(qremFit$fitted.mod))
#' summary(qremFit$fitted.mod)$coef
QREM <- function (func, linmod, dframe, qn, userwgts = NULL,..., err = 10,
                  maxit = 1000, tol = 0.001, maxInvLambda = 300) {
  ycolnum <- which(colnames(dframe) == as.character(formula(linmod)[2]))
  if (length(ycolnum) == 0)
    stop("QREM: response variable '", as.character(formula(linmod)[2]),
         "' not found in dframe.")
  y0 <- dframe[, ycolnum]
  #y0 <- model.extract(model.frame(linmod, dframe), "response")
  uwgts <- rep(1, length(y0))
  if (!is.null(userwgts)) { uwgts <- as.numeric(dframe[,userwgts]) }
  args <- list(...)
  args$weights <- uwgts
  args$formula <- linmod
  args$data <- dframe
  modelFit <- do.call(func, args )
  loglikNew <- as.numeric(logLik(modelFit))
  loglikOld <- loglikNew + 2 * err
  invLambda <- pmin(1/abs(residuals(modelFit)), maxInvLambda)
  ui <- as.vector(residuals(modelFit))
  if (inherits(modelFit, "lmerMod")) {
    modelCoefs <- list(beta = fixef(modelFit), u = ranef(modelFit))
  } else {
    modelCoefs <- list(beta = modelFit$coefficients)
  }
  it <- 0
  while ((err > tol) & ((it = it + 1) < maxit)) {
    args$weights <- invLambda*uwgts
    dframe[, ycolnum] <- y0 - (1 - 2 * qn)/invLambda
    args$data <- dframe
    modelFit <- do.call(func, args )
    if (inherits(modelFit, "lmerMod")) {
      modelFittedValues <- fitted(modelFit)
      modelCoefs <- list(beta = fixef(modelFit), u = ranef(modelFit))
    } else {
      modelFittedValues <- modelFit$fitted.values
      modelCoefs <- list(beta = modelFit$coefficients)
    }
    ui <- as.vector(y0 - modelFittedValues)
    invLambda <- pmin(1/abs(ui), maxInvLambda)
    loglikNew <- as.numeric(logLik(modelFit))
    err <- abs(loglikNew - loglikOld)
    loglikOld <- loglikNew
  }
  list(coef=modelCoefs, fitted.mod=modelFit, empq=mean(ui < 0), ui=ui,
       weights=invLambda, iter=it, err=err)
}

#' Bootstrap estimates for the standard errors of QREM coefficients.
#'
#' Runs \code{B} bootstrap replications of \code{\link{QREM}} in parallel
#' (using \code{parallel::parLapply}) and returns the full \eqn{B \times p}
#' matrix of coefficient vectors, from which standard errors can be obtained
#' as \code{apply(result, 2, sd)}.  For fixed-effects models,
#' \code{\link{bcov}} provides a faster analytic alternative based on
#' Bahadur's representation.  For mixed models, bootstrap SEs are preferred.
#'
#' @param func The fitting function (lm, lmer, gam).
#' @param linmod A formula (the linear model for fitting in the M step).
#' @param dframe0 The design matrix. A data frame containing the columns in the formula specified in linmod.
#' @param qn The selected quantile. Must be in (0,1).
#' @param n The number of rows to sample. Should equal \code{nrow(dframe0)} for a standard
#'   non-parametric bootstrap; smaller values sample a subset.
#' @param userwgts The user-provided sampling weights (optional. Default=NULL.)
#' @param ... Any arguments to be passed to func (except for the formula and weights)
#' @param sampleFrom A subset of rows in dframe0 to sample from (for mixed models). Default=NULL.
#' @param B The number of bootstrap iterations (default=100).
#' @param err The initial value for the estimation error (default=10).
#' @param maxit The maximum number of EM iterations (default=1000).
#' @param tol The error tolerance level (default=0.001).
#' @param maxInvLambda The maximum value of the weight for WLS fitting (default=300).
#' @param seedno The seed for reproducibility (default=71371).
#' @param showEst Boolean - whether to show an estimated completion time for the bootstrap. Default=FALSE.
#' @return A \eqn{B \times p} numeric matrix where each row is the fixed-effects
#'   coefficient vector from one bootstrap replicate. Column names match the
#'   coefficient names from \code{\link{QREM}}. Replications that produced a
#'   coefficient vector of unexpected length (rank-deficient resamples) are
#'   silently dropped, so the number of rows may be slightly less than \code{B}.
#' @seealso \code{\link{QREM}}, \code{\link{bcov}}
#' @importFrom parallel detectCores makeCluster parLapply clusterExport stopCluster
#' @importFrom gam gam
#' @export
#' @examples
#' \donttest{
#' # Skipped when R CMD check restricts parallel cores
#' if (!isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))) {
#'   data(simdf)
#'   estBS <- boot.QREM(lm, linmod = y ~ x * x2 + x3, dframe0 = simdf,
#'                      qn = 0.2, n = nrow(simdf), B = 50)
#'   apply(estBS, 2, sd)
#' }
#' }
boot.QREM <- function(func, linmod, dframe0, qn, n, userwgts=NULL,...,
                      sampleFrom = NULL, B = 100, err = 10,
                      maxit = 1000,tol = 0.001, maxInvLambda = 300,
                      seedno = 71371,showEst = FALSE) {
  t0 <- Sys.time()
  set.seed(seedno)
  bs_set <- sample(n)
  if (!is.null(sampleFrom)) {
    colNum <- which(colnames(dframe0) == sampleFrom)
    dframe <- dframe0[which(dframe0[, colNum] %in% bs_set), ]
  } else {
    colNum <- 0
    dframe <- dframe0[bs_set, ]
  }
  qremFit0 <- QREM(func, linmod, dframe, qn, userwgts=userwgts, ...,err=err,
                   maxit = maxit, tol = tol, maxInvLambda = maxInvLambda)
  if (B == 1) { return(qremFit0$coef$beta) }
  oneIt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  useCores <- max(detectCores() - 1, 1, na.rm = TRUE)
  if (showEst) {
    cat("One iteration ", ceiling(oneIt), "seconds\n")
    cat("Estimated completion time, using", useCores, " cores >",
        ceiling(oneIt * ceiling(B/useCores)), " seconds\n")
  }
  t0 <- Sys.time()
  n_coefs <- length(qremFit0$coef$beta)
  cl <- makeCluster(useCores)
  addArgs <- paste(list(...))
  clusterExport(cl, varlist = c("func", "linmod", "dframe0", "userwgts","qn",
                                "n","QREM", "lm", "lmer", "gam", "colNum",
                                "addArgs", "sampleFrom",
                                "fixef","ranef", "seedno", "err", "maxit",
                                "tol","maxInvLambda"), envir = environment())
  QREMpar = parLapply(cl, 1:(B - 1), function(repnum) {
    set.seed(seedno + 19 * repnum)
    bs_set <- sample(n, replace=TRUE)
    if (!is.null(sampleFrom))
      dframe <- dframe0[which(dframe0[, colNum] %in% bs_set), ]
    else dframe <- dframe0[bs_set, ]
    qremFit <- QREM(func, linmod, dframe, qn, userwgts=userwgts, ...,
                    err=err,maxit=maxit, tol=tol, maxInvLambda=maxInvLambda)
    qremFit$coef$beta
  })
  stopCluster(cl)
  if (showEst)
    cat("Actual completion time =",
        ceiling(as.numeric(difftime(Sys.time(),t0, units = "secs"))), " seconds\n")
  QREMpar <- Filter(function(x) length(x) == n_coefs, QREMpar)
  rbind(qremFit0$coef$beta, matrix(unlist(QREMpar), ncol = n_coefs,
                                   byrow = TRUE))
}


#' Covariance matrix estimation
#'
#' Estimates the covariance matrix of the fixed-effects coefficients using
#' Bahadur's representation:
#' \deqn{Cov(\hat{\beta}) = \frac{\tau(1-\tau)}{f(0)^2} (X'X)^{-1}}
#' where \eqn{f(0)} is the density of the residuals at zero, estimated via a
#' kernel density estimator using Silverman's bandwidth
#' (\eqn{0.9 \cdot \min(s, IQR/1.34) \cdot n^{-0.2}}).
#' In the mixed model case, the random effects are treated as fixed (BLUPs),
#' and the same formula is applied; only the fixed-effects sub-matrix is returned.
#' @param qremFit A fitted model object returned by \code{\link{QREM}}.
#' @param linmod A formula (the linear model for fitting in the M step).
#' @param dframe A data frame containing the columns in the formula.
#' @param qn The selected quantile. Must be in (0,1).
#' @param userwgts A column name or index in \code{dframe} for user-provided
#'   sampling weights (optional. Default=NULL). Sampling weights are applied in
#'   the fixed-effects case only; ignored for \code{lmerMod} fits.
#' @return A covariance matrix for the fixed-effects coefficients. For mixed
#'   models, only the fixed-effects sub-matrix (rows/columns corresponding to
#'   \code{X}) is returned, not the full augmented matrix including random effects.
#' @seealso \code{\link{QREM}}, \code{\link{boot.QREM}}
#' @importFrom KernSmooth bkde
#' @importFrom stats IQR fitted formula logLik model.extract model.frame model.matrix qqplot quantile residuals sd as.formula lm
#' @importFrom lme4 getME
#' @importFrom corpcor make.positive.definite
#' @export
#' @examples
#' data(simdf)
#' qremFit <- QREM(lm, linmod = y ~ x * x2 + x3, dframe = simdf, qn = 0.2)
#' covmat <- bcov(qremFit, linmod = y ~ x * x2 + x3, dframe = simdf, qn = 0.2)
bcov <- function (qremFit, linmod, dframe, qn, userwgts=NULL) {
  ui <- qremFit$ui
  empq <- quantile(ui, qn)
  ui <- ui - empq
  if (!is.null(userwgts)) { ui <- ui*sqrt(userwgts) }
  # estimate f(0) with a kernel density estimator
  bw <- (length(ui))^(-0.2) * 0.9 * min(sd(ui), IQR(ui)/1.34)
  kdest <- bkde(ui, bandwidth = bw)
  negIdx <- which(kdest$x < 0)
  if (length(negIdx) == 0 || negIdx[length(negIdx)] >= length(kdest$x)) {
    warning("bcov: residuals have no mass near zero; f(0) estimate is unreliable.")
    return(NULL)
  }
  largestNeg <- negIdx[length(negIdx)]  # index in full kdest$x of the last negative x
  xs <- c(kdest$x[largestNeg + 1], -kdest$x[largestNeg])
  ys <- c(kdest$y[largestNeg], kdest$y[largestNeg + 1])
  dFp0 <- sum(xs * ys)/sum(xs)
  if (dFp0 < .Machine$double.eps) {
    warning("bcov: estimated density at zero is near zero; covariance estimates will be unreliable.")
  }
  if (inherits(qremFit$fitted.mod, "lmerMod")) {
    XX <- getME(qremFit$fitted.mod,"X")
    ZZ <- getME(qremFit$fitted.mod,"Z")
    DM <- if (ncol(ZZ) > 1) cbind(XX, ZZ[,-ncol(ZZ)]) else XX
  } else {
    DM <- XX <- model.matrix(linmod, dframe)
    if (!is.null(userwgts)) { # in the fixed effect model only
      DM <- Diagonal(length(userwgts),userwgts)%*%DM
    }
  }
  invMat <- solve(make.positive.definite(t(DM) %*% DM))
  as.matrix(invMat * qn * (1 - qn)/(dFp0^2))[1:ncol(XX), 1:ncol(XX)]
}


#' Quantile regression diagnostics based on the residuals, u_i.
#'
#' For a given predictor, create a qq-plot if continuous, or a spinogram if categorical.
#'
#' @param X A predictor included in the regression.
#' @param varname The predictor's name (for plotting).
#' @param u_i The QR residuals.
#' @param qn The quantile used in the regression.
#' @param plot.it Boolean, if TRUE, will show the histogram of the residuals and the fitted kernel density estimate (default=TRUE).
#' @param filename The pdf file to save the plot. Default is NULL (print to the screen.)
#' @return For a continuous predictor: a list with elements \code{x}, \code{y}
#'   (the QQ-plot quantile vectors from \code{\link[stats]{qqplot}}) and
#'   \code{dev} (the marginal quantile deviance,
#'   \eqn{-\sum u_i (\tau - \mathbf{1}_{u_i > 0})}).
#'   For a factor predictor: a named list with one numeric entry per factor
#'   level (the empirical proportion of observations below the regression line
#'   for that level) plus \code{dev}.
#'   Returns \code{NULL} with a \code{warning} if residuals are constant, all
#'   fall on one side of zero, or \code{X} is neither numeric nor factor.
#' @importFrom graphics abline axis grid hist lines plot rect text
#' @export
#' @examples
#' data(simdf)
#' qremFit <- QREM(lm, linmod = y ~ x * x2 + x3, dframe = simdf, qn = 0.2)
#' qrdg <- QRdiagnostics(simdf$x, "x", qremFit$ui, 0.2, plot.it = FALSE)
#' qrdg <- QRdiagnostics(simdf$x3, "x3", qremFit$ui, 0.2, plot.it = FALSE)
QRdiagnostics <- function(X, varname, u_i, qn,  plot.it=TRUE, filename=NULL) {
  n <- length(u_i)
  empq <- quantile(u_i,qn)
  u_i <- u_i-empq
  if (sd(u_i) == 0) {
    warning("QRdiagnostics: all residuals are identical; cannot estimate density.")
    return(NULL)
  }
  kde <- bkde(u_i)
  # calculate the deviance:
  mod.deviance <- -sum(u_i*(qn-(u_i >0)))
  # show the histogram of the residuals and the fitted kernel density estimate:
  if(plot.it && is.null(filename)) {
    hist(u_i,breaks=ifelse(n > 500, 40, 10),freq=FALSE)
    lines(kde, col=2)
    abline(v=0,col=4,lwd=2)
  }
  # for a given predictor, create a qq-plot (continuous) or a spinogram (categorical)
  posidx <- which(u_i > 0)
  if (length(posidx) == 0 || length(posidx) == n) {
    if (plot.it && !is.null(filename)) dev.off()
    warning("QRdiagnostics: all residuals are on the same side of zero; qq-plot is not meaningful.")
    return(NULL)
  }
  if (plot.it && ! is.null(filename))
    pdf(filename, width=5, height=5)
  if (!is.numeric(X) && !is.factor(X)) {
    if (plot.it && !is.null(filename)) dev.off()
    warning("QRdiagnostics: X must be numeric or factor; no plot produced.")
    return(NULL)
  }
  if (is.numeric(X)) {
    qqp <- qqplot(X[posidx] , X[-posidx], xlim=c(min(X),max(X)),
                  ylim=c(min(X),max(X)),pch=19,cex=0.4,col="blue",
                  xlab = "Q(above)", ylab="Q(below)",
                  main=paste(varname,", q=",qn),  plot.it = plot.it)
    if(plot.it) {
      abline(0,1,col=2, lwd=2)
      grid()
    }
    if (plot.it && ! is.null(filename))
      dev.off()
    qqp$dev <- mod.deviance
    return(qqp)
  }
  if (is.factor(X)) {
    qqlvl <- as.list(colSums(table(u_i[which(u_i < 0)], X[which(u_i < 0)])) / table(X))
    if(plot.it) {
      plot(c(0.25,length(levels(X)))+0.5,c(0,1.2), col=0, axes=F,
           main = paste(varname,", q=",qn), xlab="", ylab="Prob(u < 0)")
      text(1:length(levels(X)), 1.1, levels(X))
      for (i in 1:length(levels(X))) {
        rect(i-0.25,0,i+0.25,qqlvl[i], col = "grey33", border="white")
        rect(i-0.25,qqlvl[i],i+0.25,1, col = "grey77", border="white")
        lines(c(i-0.25, i+0.25), c(qn, qn), col=2, lwd=3)
      }
      abline(h=seq(0.1,0.9, by=0.1), col="white", lwd=0.5, lty=3)
    }
    if (plot.it && ! is.null(filename))
      dev.off()
    qqlvl$dev <- mod.deviance
    return(qqlvl)
  }
}



#' A `flat QQ-plot' for a fitted quantile regression model.
#'
#' Showing multiple ('flat') QQ plots for different quantiles using a heatmap.
#'
#' @param dat The data matrix.
#' @param cnum The selected column number.
#' @param vname The selected variable name (if cnum is not specified).
#' @param qrfits A list of QR-fitted models.
#' @param qns The quantiles.
#' @param maxm The maximum number of segments in the partision of the selected variable (default=30)
#' @param sdevs Used to determine the ranges for critical values under the null model and thus the range of colors in the heatmap (default= 4).
#' @param filename The pdf file to save the plot. Default is NULL (print to the screen.)
#' @param plot.it A logical variable used to determine whether to show the diagnostic plot (TRUE).
#' @return A list of length \code{length(qns)}, where each element is an
#'   \code{htest} object from \code{\link[stats]{prop.test}} testing whether
#'   the proportion of observations below the regression line equals
#'   \code{qns[i]} in each segment of the predictor. Returns \code{FALSE}
#'   (with a \code{warning}) if neither \code{cnum} nor \code{vname} is
#'   specified, or if the sample size is too small.
#' @seealso \code{\link{QREM}}, \code{\link{QRdiagnostics}}
#' @importFrom grDevices topo.colors pdf dev.off
#' @importFrom stats prop.test qnorm
#' @export
#' @examples
#' \donttest{
#' data(simdf)
#' qns <- seq(0.1, 0.9, by = 0.1)
#' qrfits <- lapply(qns, function(q)
#'   QREM(lm, linmod = y ~ x * x2 + x3, dframe = simdf, qn = q))
#' pvals <- flatQQplot(dat=simdf, cnum = 2, qrfits=qrfits, qns=qns, maxm = 20, plot.it = TRUE)
#' pvals <- flatQQplot(dat=simdf, cnum = 3, qrfits=qrfits, qns=qns, maxm = 20, plot.it = TRUE)
#' pvals <- flatQQplot(dat=simdf, cnum = 4, qrfits=qrfits, qns=qns, maxm = 20, plot.it = TRUE)
#' }
flatQQplot <- function(dat, cnum=NULL, vname=NULL, qrfits, qns,
                          maxm=30, sdevs = 4, filename = NULL, plot.it=TRUE) {
  if (is.null(cnum) & is.null(vname)) {
    warning("Must specify either a valid column number or a column name.")
    return(FALSE)
  }
  if (!is.null(vname)) {
    cnum <- which(colnames(dat) == vname)
  }
  x <- dat[,cnum]

  N <- nrow(dat)
  k <- length(qns)
  if (is.factor(dat[,cnum])) {
    x <- droplevels(x)
    m <- length(levels(x))
    endpts <- 1:length(levels(x))+0.5
    xcoord <- 1:length(levels(x))
    xcoord <- c(xcoord, max(xcoord)+1)
  }
  if (is.numeric(dat[,cnum])) {
    m <- min(floor(N/(10*k)), maxm)
    if (m <= 1) {
      warning("Sample size is not large enough to perform this function for ", k, " quantiles")
      return(FALSE)
    }
    endpts <- c(quantile(x, probs = (1:(m-1))/m), max(x))
    names(endpts) <- c()
    xcoord <- unique(c(min(x), endpts))
    m <- length(xcoord) - 1
  }
  obscnt <- sqcols <- matrix(0, ncol=m, nrow=k)
  pvals <- list()
  cutoffs <- qnorm(c(0, kronecker(10^seq(-sdevs, -1), c(0.25, 0.5, 1))))
  cutoffs <- c(cutoffs, rev(-cutoffs))
  L <- length(cutoffs)
  cols <- topo.colors(L)
  for (i in 1:k) {
    qremFit <-  qrfits[[i]]
    belowQR <- (qremFit$ui < 0)
    if (is.numeric(dat[,cnum])) {
      Gij <- cut(x, breaks = xcoord, include.lowest = TRUE)
      gijtab <- table(belowQR, Gij)
      if (nrow(gijtab) == 1) {
        if(rownames(gijtab) == "FALSE") {
          obscnt[i,] <- rep(0, ncol(gijtab))
        } else {
          obscnt[i,] <- gijtab[1,]
        }
      } else {
        obscnt[i,] <- gijtab[2,]
      }
    }
    if (is.factor(dat[,cnum])) {
      Gij <- x
      gijtab <- table(belowQR, dat[,cnum])
      if (nrow(gijtab) == 1) {
        if (rownames(gijtab) == "FALSE") {
          obscnt[i,] <- rep(0, ncol(gijtab))
        } else {
          obscnt[i,] <- gijtab[1,]
        }
      } else {
        obscnt[i,] <- gijtab[2,]
      }
    }
    pvals[[i]] <- prop.test(obscnt[i,], pmax(table(Gij), 0.1), rep(qns[i], m))
    sqcols[i,] <- as.numeric(cut((obscnt[i,]/table(Gij)-qns[i])/sqrt(qns[i]*(1-qns[i])/table(Gij)),breaks = cutoffs))
  }
  mdist <- 0.5*max(abs(diff(xcoord)))
  qdist <- ifelse(length(qns) > 1, 0.5*max(abs(diff(qns))),1)
  if (!is.null(filename))
    pdf(filename, width = 5, height = 5)
  if (plot.it) {
    plot( c(min(xcoord) - 0.5 * mdist, max(xcoord) + 5 * mdist),
          c(0 - qdist, 1),  col = 0, axes = F,
          main = colnames(dat)[cnum], ylab = "quantile", xlab = "x")
    if (is.numeric(dat[,cnum]))
      axis(1, at = sprintf("%3.1f",xcoord))
    if (is.factor(dat[,cnum]))
      text(xcoord[1:m]+mdist,rep(0,m), levels(x))
    text(min(xcoord),qns,  qns, cex = 0.7, pos=2)
    #  Chi <- (obscnt-expcnt)/sqrt(expcnt)
    for (j in 1:ncol(obscnt)) {
      for (i in 1:nrow(obscnt)) {
        rect( xcoord[j], qns[i] - 0.95*qdist,
              xcoord[j + 1], qns[i] + 0.95*qdist,
              col = if (is.na(sqcols[i,j])) cols[L] else cols[sqcols[i,j]],
              border = "white")
      }
    }
    # legend:
    #rngx <- max(qns) - min(qns)
    rngx <- 1
    for (j in 1:L) {
      rect(max(xcoord) + 0.25 * mdist, (j - 1)*rngx/L,
           max(xcoord) + mdist, j*rngx/L, col = cols[j], border = "white")
      text(max(xcoord) + 2*mdist, (j-0.5)*rngx/L,
           format(cutoffs[j], digits=3), cex=0.6, col="grey33", font = 2)
    }
  }
  if (!is.null(filename))
    dev.off()
  return(pvals)
}



#' Variable selection for quantile regression.
#'
#' Use the SEMMS package to perform variable selection. Iteratively alternate between QREM and fitSEMMS.
#'
#' @param inputData A data frame or a file generated for a SEMMS analysis. See the SEMMS package for details. If a data frame, it will be saved in a tempfile.
#' @param ycol The number of the column in the input file which should be used as the response.
#' @param Zcols The columns in the input file which contain the putative variables.
#' @param Xcols The columns in the input file which contain the fixed effects in the model (default is none, c()).
#' @param qn The selected quantile. Must be in (0,1).
#' @param nn The initial value for the number of non-null variables in SEMMS. Default is 5.
#' @param nnset Optional: instead of an initial number of candidates, can specify the column numbers in the Z matrix for the first iteration. Default is null.
#' @param maxRep The maximum number of iterations between QREM and fitSEMMS. Default=40.
#' @param initWithEdgeFinder Determines whether to use the edgefinder package to find highly correlated pairs of predictors (default=FALSE).
#' @param mincor To be passed to the fitSEMMS function (the minimum correlation coefficient between pairs of putative variable, over which they are considered highly correlated). Default is 0.75.
#' @return A list with two elements:
#' \itemize{
#'   \item \code{fittedSEMMS}: The SEMMS model object from the final iteration
#'     (output of \code{fitSEMMS}), containing the selected predictor indices
#'     in \code{$gam.out$nn}.
#'   \item \code{fittedQREM}: The QREM fit from the final iteration
#'     (output of \code{\link{QREM}}), or \code{NULL} if SEMMS selected no
#'     predictors.
#' }
#' @seealso \code{\link{QREM}}, \code{\link{bcov}}
#' @importFrom SEMMS fitSEMMS readInputFile
#' @importFrom edgefinder edgefinder graphComponents
#' @importFrom Matrix Diagonal
#' @export
#' @examples
#' \donttest{
#' data(simLargeP)
#' qn <- 0.25
#' res <- QREM_vs(simLargeP, 1, 2:51, qn=qn)
#' dfsemms <- simLargeP[,c(1, 1+res$fittedSEMMS$gam.out$nn)]
#' qremFit <- QREM(lm, y~., dfsemms, qn=qn)
#' ests <- rbind(qremFit$coef$beta,
#'          sqrt(diag(bcov(qremFit, linmod = y ~ ., dframe = dfsemms, qn = qn))))
#' rownames(ests) <- c("Estimate","s.d")
#' }
QREM_vs <- function(inputData, ycol, Zcols, Xcols=c(), qn, nn=5, nnset=NULL, maxRep=40,initWithEdgeFinder=FALSE, mincor = 0.75) {
  if (!is.data.frame(inputData) && !is.character(inputData)) {
    stop("Invalid input type: ", class(inputData), "\n")
  }
  if (is.data.frame(inputData)) {
    filename <- tempfile(pattern = "forsemms_",fileext = ".RData")
    save(inputData, file=filename)
  } else {
    filename <- inputData
  }
  if(!file.exists(filename)) {
    stop("File doesn't exist:", filename,"\n")
  }
  dataYXZ <- readInputFile(filename, ycol=ycol, Xcols = Xcols, Zcols=Zcols)
  y0 <- scale(dataYXZ$Y)
  if(initWithEdgeFinder) {
    M <- t(cbind(y0, dataYXZ$Z))
    effit <- edgefinder(M, BHthr = max(1e-6, 1/choose(nrow(M),2)), LOvals = 100)
    subgr <- graphComponents(effit$AdjMat[-1,-1], minCtr = 2)
    Zcols <- sort(c(which(subgr$clustNo == 0), which(subgr$iscenter ==1)))
    if (!is.null(nnset)) {
      Zcols <- sort(union(nnset, Zcols))
      nnset <- which(Zcols %in% nnset)
    } else {
      if(length(which(effit$AdjMat[1,setdiff(Zcols,1)] != 0) > 0))
        nnset <- which(effit$AdjMat[1,setdiff(Zcols,1)] != 0)
    }
    dat0 <- dataYXZ
    dataYXZ$Z <- dat0$Z[,Zcols]
    dataYXZ$K <- length(Zcols)
    dataYXZ$originalZnames <- dat0$originalZnames[Zcols]
    dataYXZ$colnamesZ <- dat0$colnamesZ[Zcols]
  }
  if (is.null(nnset)) { # get the initial set of predictors, if not provided
    zval <- rep(0, dataYXZ$K)
    rnd <- sample(dataYXZ$K, replace=FALSE)
    m <- 5
    for (i in 1:ceiling(dataYXZ$K/m)) {
      idx <- ((i-1)*m+1) : min(i*m, dataYXZ$K)
      preds <- paste0(colnames(dataYXZ$Z)[rnd[idx]], collapse = " + ")
      linmod <- as.formula(paste("Y ~", preds))
      dframetmp <- data.frame(cbind(dataYXZ$Y, dataYXZ$Z[,rnd[idx]]))
      colnames(dframetmp) <- c("Y",colnames(dataYXZ$Z)[rnd[idx]])
      qremFit <- QREM(lm, linmod, dframetmp, qn, maxInvLambda = 1000)
      se <- sqrt(pmax(diag(bcov(qremFit,linmod,dframetmp,qn)), 0))[-1]
      zval[rnd[idx]] <- ifelse(se > 0, qremFit$coef$beta[-1]/se, 0)
    }
    nnset <- order(abs(zval),decreasing = TRUE)[1:nn]
  }
  prev_ll <- 0
  fittedVSnew <- NULL
  for (repno in 1:maxRep) {
    # create a subset of the selected columns and run QREM
    preds <- paste(colnames(dataYXZ$Z)[nnset], collapse = "+")
    linmod <- as.formula(paste("Y ~", preds))
    dframetmp <- data.frame(cbind(dataYXZ$Y, dataYXZ$Z[,nnset]))
    colnames(dframetmp) <- c("Y", colnames(dataYXZ$Z)[nnset])
    qremFit <- QREM(lm, linmod, dframetmp, qn, maxInvLambda = 1000)
    new_ll <- sum(qremFit$ui*(qn - as.numeric(qremFit$ui < 0)))
    if (abs(prev_ll - new_ll) < 1e-3) {
      if(initWithEdgeFinder) {
        A <- effit$AdjMat
        if (length(fittedVSnew$gam.out$nn) > 0){
          fittedVSnew$gam.out$nn <- Zcols[fittedVSnew$gam.out$nn]
          A[1, fittedVSnew$gam.out$nn+1] <- A[fittedVSnew$gam.out$nn+1, 1] <- 1
        }
        fittedVSnew$gam.out$lockedOut <- rep(0, dat0$K)
        lout <- setdiff(which((A + A%*%A)[1,] > 0), 1)
        fittedVSnew$gam.out$lockedOut[setdiff(lout-1, fittedVSnew$gam.out$nn)] <- 1
        fittedVSnew$inits$discard <- rep(-1, dat0$K)
        fittedVSnew$gam.out$A <- effit$AdjMat
      }
      return(list(fittedSEMMS=fittedVSnew, fittedQREM=qremFit))
    }
    prev_ll <- new_ll
    # apply the weights found by QREM and run SEMMS
    dataYXZtmp <- dataYXZ
    dataYXZtmp$Y <- (dataYXZ$Y - (1-2*qn)/qremFit$weights)
    fittedVSnew <- fitSEMMS(dataYXZtmp, distribution = 'N', mincor=mincor, rnd=F,
                            nnset=nnset, minchange = 1, maxst = 20,
                            initWithEdgeFinder=FALSE)
    if(initWithEdgeFinder) {
      A <- effit$AdjMat
      if (length(fittedVSnew$gam.out$nn) > 0) {
        fittedVSnew$gam.out$nn <- Zcols[fittedVSnew$gam.out$nn]
        A[1, fittedVSnew$gam.out$nn+1] <- A[fittedVSnew$gam.out$nn+1, 1] <- 1
      }
      fittedVSnew$gam.out$lockedOut <- rep(0, dat0$K)
      lout <- setdiff(which((A + A%*%A)[1,] > 0), 1)
      fittedVSnew$gam.out$lockedOut[setdiff(lout-1, fittedVSnew$gam.out$nn)] <- 1
      fittedVSnew$inits$discard <- rep(-1, dat0$K)
      fittedVSnew$gam.out$A <- effit$AdjMat
    }
    if (length(fittedVSnew$gam.out$nn) == 0)
      return(list(fittedSEMMS=fittedVSnew, fittedQREM=NULL))
    nnset <- fittedVSnew$gam.out$nn
  }
  return(list(fittedSEMMS=fittedVSnew, fittedQREM=qremFit))
}
