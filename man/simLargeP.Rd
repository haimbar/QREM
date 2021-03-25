\name{simLargeP}
\alias{simLargeP}
\docType{data}
\title{
simLargeP
}
\description{
A simulated data set for the QREM package with 500 predictors, of which only 5 are associated with the outcome. The standard deviation of the random errors is a linear function of the the first predictor.
}
\usage{data("simLargeP")}
\format{
  A data frame with 200 simulated observations with the 500 variables, of which the first 5 are related to the mean of y. The coefficients, including the intercept are c(1, -3,2,2,-1,-2). The errors are generated using a normal distribution with standard deviation equal to 0.1+X1.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
# See in the description of the QREM_vs function.
}
\keyword{datasets}
