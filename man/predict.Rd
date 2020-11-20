\encoding{UTF-8}
\name{predict}
\alias{predict}
\alias{predict.default}
\title{Forecasting with Singular Spectrum Trendline}
\description{
  Computes a forcasted trendline for time series data using
  singular spectrum analysis. 
}
\usage{
predict(fitted.model, p = 1)
}
\arguments{
  \item{fitted.model}{estimated model using the functions in the package.}
  \item{p}{the horizon to produce forecasts.}
  }
\value{
  \item{forecast}{Matrix containing in columns the dimensions and in rows the forecasts.}
  \item{a}{Parameters corresponding to the linear recurrence formula.}
}
\details{
  predict is a wrapper function for predictions from the results of various singular spectrum model fitting functions on the package. The function invokes particular methods which depend on the class of the first argument (i.e. \link{sst},\link{msst}, \link{msstc}, \link{isst}, \link{misst}, etc) .
}
\seealso{
  See \code{\link{sst},\link{msst}, \link{msstc}, \link{isst}, \link{misst}} 
  for a version of   different models.
} 
\references{
  de Carvalho, M. and Martos, G. (2020). Brexit: Tracking and
  disentangling the sentiment towards leaving the EU. \emph{International Journal
  of Forecasting}, \bold{36}, 1128--1137.
  de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
}
\examples{
## SIMULATED DATA EXAMPLE
set.seed(1)
N <- 500
t <- seq(.01, 5, length = N)
Y <- cbind(t^3 - 9 * t^2 + 23 * t + rnorm(N, 0, 1), 
           10 * sin(3 * t) / t + rnorm(N, 0, 1))
y <- mtsframe(date = t, Y)
fit <- msst(y)

# Forecasting: 
predict(fit, p = 5)$forecast
}
\keyword{datasets}
\author{Gabriel Martos and Miguel de Carvalho}
