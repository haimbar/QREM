#' Fitting a quantile regression model via the EM algorithm.
#'
#' @param func The fitting function (lm, lmer, or gam).
#' @param linmod A formula (the linear model for fitting in the M step).
#' @param dframe A data frame containing the columns in the formula.
#' @param qn The selected quantile. Must be in (0,1).
#' @param userwgts The user-provided sampling weights (optional. Default=NULL.)
#' @param ... Any arguments to be passed to func (except for the formula and weights). Note that gam requires the family of the error distribution.
#' @param err The initial value for the estimation error (default=10). Must be greater than tol (below).
#' @param maxit The maximum number of EM iterations (default=1000).
#' @param tol The error tolerance level (default=0.001).
#' @param maxInvLambda The maximum value of the weight for WLS fitting (default=300).
#' @return A list containing the following
#' \itemize{
#'   \item coef The estimated regression coefficients (a list).
#'   \item fitted.mod The output from lm(), lmer(), or gam() in the last EM iteration.
#'   \item empq The percentage of points below the regression line (the empirical quantile).
#'   \item ui The quantile regression residuals.
#'   \item weights The weights used in the WLS solution.
#'   \item iter The number of EM iterations.
#'   \item err The final EM error (convergence criterion).
#' }
#' @importFrom lme4 lmer fixef ranef
#' @export
#' @examples
#' data(simdf)
#' qremFit <-  QREM(lm,linmod=y~x*x2 +x3, df=simdf, qn=0.2)
#' summary(aov(qremFit$fitted.mod))
#' summary(qremFit$fitted.mod)$coef
QREM <- function (func, linmod, dframe, qn, userwgts = NULL,..., err = 10,
                  maxit = 1000, tol = 0.001, maxInvLambda = 300) {
  ycolnum <- which(colnames(dframe) == as.character(formula(linmod)[2]))
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
  it <- 0
  while ((err > tol) & ((it = it + 1) < maxit)) {
    args$weights <- invLambda*uwgts
    dframe[, ycolnum] <- y0 - (1 - 2 * qn)/invLambda
    args$data <- dframe
    modelFit <- do.call(func, args )
    if (class(modelFit) == "lmerMod") {
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

#' Bootstrap estimates for the standard errors of the coefficients in a quantile regression model.
#'
#' In the fixed effects case, the bcov function provides less variable,
#'  and faster estimates through the asymptotic covariance
#'   (Bahadur's representation). For mixed models bcov may also be used - it provides
#'   good coverage probability in simulations (using the BLUPs for the random
#'   effects)
#'
#' @param func The fitting function (lm, lmer, gam).
#' @param linmod A formula (the linear model for fitting in the M step).
#' @param dframe0 The design matrix. A data frame containing the columns in the formula specified in linmod.
#' @param qn The selected quantile. Must be in (0,1).
#' @param n The number of samples to be used in the bootstrap.
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
#' @return A matrix of the QR coefficients (B rows).
#' @importFrom parallel detectCores makeCluster parLapply clusterExport stopCluster
#' @importFrom gam gam
#' @export
#' @examples
#' \donttest{
#' #data(simdf)
#' #qremFit <-  QREM(lm,linmod=y~x*x2 +x3, df=simdf, qn=0.2)
#' #estBS <- boot.QREM(lm, linmod=y~x*x2 +x3, df = simdf, qn=0.2,
#' #    n=nrow(simdf), B=50)
#' #apply(estBS,2,sd)
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
  useCores <- detectCores() - 1
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
                                "addArgs",
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
  rbind(qremFit0$coef$beta, matrix(unlist(QREMpar), ncol = n_coefs,
                                   byrow = TRUE))
}


#' Covariance matrix estimation
#'
#' Covariance matrix estimation for the fixed effects case using Bahadur's
#' representation. In mixed model, we use the BLUPs (as if they were estimated
#' as fixed effects) and the same formula is used.
#' @param qremFit A fitted model object, returned by QREM.
#' @param linmod A formula (the linear model for fitting in the M step).
#' @param dframe A data frame containing the columns in the formula.
#' @param qn The selected quantile. Must be in (0,1).
#' @param userwgts The user-provided sampling weights (optional. Default=NULL.)
#' @return A covariance matrix for the fixed-effects model.
#' @importFrom KernSmooth bkde
#' @importFrom stats IQR fitted formula logLik model.extract model.frame model.matrix qqplot quantile residuals sd as.formula lm
#' @importFrom lme4 getME
#' @importFrom corpcor make.positive.definite
#' @export
#' @examples
#' data(simdf)
#' qremFit <-  QREM(lm,linmod=y~x*x2 +x3, df=simdf, qn=0.2)
#' covmat <- bcov(qremFit,linmod=y~x*x2 +x3, df=simdf, qn=0.2)
bcov <- function (qremFit, linmod, dframe, qn, userwgts=NULL) {
  ui <- qremFit$ui
  empq <- quantile(ui, qn)
  ui <- ui - empq
  if (!is.null(userwgts)) { ui <- ui*sqrt(userwgts) }
  # estimate f(0) with a kernel density estimator
  bw <- (length(ui))^(-0.2) * 0.9 * min(sd(ui), IQR(ui)/1.34)
  kdest <- bkde(ui, bandwidth = bw)
  largestNeg <- which.max(kdest$x[kdest$x < 0])
  xs <- c(kdest$x[largestNeg + 1], -kdest$x[largestNeg])
  ys <- c(kdest$y[largestNeg], kdest$y[largestNeg + 1])
  dFp0 <- sum(xs * ys)/sum(xs)
  if (class(qremFit$fitted.mod) == "lmerMod") {
    XX <- getME(qremFit$fitted.mod,"X")
    ZZ <- getME(qremFit$fitted.mod,"Z")
    DM <- cbind(XX, ZZ[,-ncol(ZZ)])
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
#' @return A list, as follows, plus the marginal deviance:
#' \itemize{
#'   \item For a continuous predictor, the list is called qqp and contains the output from qqplot().
#'   \item For a categorical variable, the list is called qqlvl, and it contains the empirical percentages of points below the regression line, for each level.
#' }
#' @importFrom graphics abline axis grid hist lines plot rect text
#' @export
#' @examples
#' data(simdf)
#' qremFit <-  QREM(lm,linmod=y~x*x2 +x3, df=simdf, qn=0.2)
#' qrdg <- QRdiagnostics(simdf$x, "x",qremFit$ui, 0.2,  plot.it = TRUE)
#' qrdg <- QRdiagnostics(simdf$x3, "x3",qremFit$ui, 0.2,  plot.it = TRUE)
QRdiagnostics <- function(X, varname, u_i, qn,  plot.it=TRUE, filename=NULL) {
  n <- length(u_i)
  empq <- quantile(u_i,qn)
  u_i <- u_i-empq
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
  if (plot.it && ! is.null(filename))
    pdf(filename, width=5, height=5)
  if(class(X) == "numeric") {
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
  if(class(X) == "factor") {
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
#' @importFrom grDevices topo.colors pdf dev.off
#' @importFrom stats prop.test qnorm
#' @export
#' @examples
#' \donttest{
#' data(simdf)
#' qns <- seq(0.1,0.9,by=0.1)
#' qrfits <- list()
#' for (i in 1:length(qns)) {
#'  qrfits[[i]] <- QREM(lm,linmod=y~x*x2 +x3, df=simdf, qn=qns[i])
#' }
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
  if (class(dat[,cnum]) == "factor") {
    x <- droplevels(x)
    m <- length(levels(x))
    endpts <- 1:length(levels(x))+0.5
    xcoord <- 1:length(levels(x))
    xcoord <- c(xcoord, max(xcoord)+1)
  }
  if (class(dat[,cnum]) %in%  c("integer","numeric")) {
    m <- min(floor(N/(10*k)), maxm)
    if (m == 1) {
      warnings("Sample size is not large enough to perform this function for ", k, "quantiles")
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
    if (class(dat[,cnum]) %in% c("integer","numeric")) {
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
    if (class(dat[,cnum]) == "factor") {
      Gij <- x
      obscnt[i,] <- table(belowQR, dat[,cnum])[2,]
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
    if (class(dat[,cnum]) %in% c("integer","numeric"))
      axis(1, at = sprintf("%3.1f",xcoord))
    if (class(dat[,cnum]) == "factor")
      text(xcoord[1:m]+mdist,rep(0,m), levels(x))
    text(min(xcoord),qns,  qns, cex = 0.7, pos=2)
    #  Chi <- (obscnt-expcnt)/sqrt(expcnt)
    for (j in 1:ncol(obscnt)) {
      for (i in 1:nrow(obscnt)) {
        rect( xcoord[j], qns[i] - 0.95*qdist,
              xcoord[j + 1], qns[i] + 0.95*qdist,
              col = cols[sqcols[i,j]],border = "white")
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
#'          sqrt(diag(bcov(qremFit,linmod=y~., df=dfsemms, qn=qn))))
#' rownames(ests) <- c("Estimate","s.d")
#' }
QREM_vs <- function(inputData, ycol, Zcols, Xcols=c(), qn, nn=5, nnset=NULL, maxRep=40,initWithEdgeFinder=FALSE, mincor = 0.75) {
  if(!(class(inputData) %in% c("data.frame","character"))) {
    stop("Invalid input type:",class(inputData), "\n")
  }
  if (class(inputData) == "data.frame") {
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
    for (i in 1:(dataYXZ$K/m)) {
      idx <- ((i-1)*m+1) : (i*m)
      preds <- paste0(colnames(dataYXZ$Z)[rnd[idx]], collapse = " + ")
      linmod <- as.formula(paste("Y ~", preds))
      dframetmp <- data.frame(cbind(dataYXZ$Y, dataYXZ$Z[,rnd[idx]]))
      colnames(dframetmp) <- c("Y",colnames(dataYXZ$Z)[rnd[idx]])
      qremFit <- QREM(lm, linmod, dframetmp, qn, maxInvLambda = 1000)
      zval[rnd[idx]] <- qremFit$coef$beta[-1]/sqrt(diag(bcov(qremFit,linmod,dframetmp,qn)))[-1]
    }
    nnset <- order(abs(zval),decreasing = TRUE)[1:nn]
  }
  prev_ll <- 0
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
