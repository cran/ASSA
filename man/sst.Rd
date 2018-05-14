\encoding{UTF-8}
\name{sst}
\alias{sst}
\alias{sst.default}
\title{Singular Spectrum Trendline}
\description{
  Computes a trendline for univariate time series data using
  singular spectrum analysis. 
}
\usage{
sst(y, l = "automatic", m = "automatic")
}
\arguments{
  \item{y}{mtsframe format data containing data. While \code{y}
    can include several time series, a more appropriate method for
    multivariate time series is \code{\link{msst}}.}
  \item{l}{window length; the string \code{"automatic"} automatic sets
    the default option \cr \code{l = ceiling(y$n + 1) / 2}.}
  \item{m}{number of leading eigentriples; the string \code{"automatic"}
    yields a vector containing the number of components in each
    dimension to be used in the trendline estimation. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default; see details.}}
\value{
  \item{trendline}{mtsframe object with trendline estimation from
    targeted grouping based on a cumulative periodogram criterion (or
    according to the number of components specified in vector
    \code{m}).}
  \item{residuals}{mtsframe object with the residuals from targeted
    grouping based on a cumulative periodogram criterion (or according
    to the number of components specified in vector \code{m}).}
  \item{erc}{list with elementary reconstructed components.}
  \item{eigen.val}{vector with the singular values of the
    trajectory matrix.}
  \item{l}{window length.}
  \item{selected.components}{vector with components selected on each
    dimension.}
  \item{selection.criteria}{a vector indicating if the hypothesis of
    white noise residual is rejected along the dimensions (0: not
    rejected, 1: rejected).}
  \item{rank}{rank of the trajectory matrix.}
}
\details{
  Singular spectrum analysis decompose time series data (\code{y})
  into principal components, and a cumulative periodogram-based
  criterion learn about elementary reconstructed components
  (\code{erc}) that contribute to the signal. The trendline results
  from adding principal components selected by a cumulative
  periodogram-based criteria; see de Carvalho and Martos (2018,
  Section 4.1). The \code{plot} method yields the resulting trendlines
  along with the data; \code{options} for the plot method are give by
  a list including the strings \code{"trendlines"},
  \code{"components"}, \code{"cpgrams"}, and \code{"screeplots"},
  along with a set of values (\code{ncomp}) indicating the components
  on which these diagnostics are to be depicted (e.g. \code{plot(fit,
  options = list(type = "components", ncomp = 1:3))}.
}
\seealso{
  See \code{\link{msst}} for a version of the routine for multivariate
  time series, and see \code{\link{msstc}} for a version of the routine
  for multivariate time series of compositional data.
} 
\references{
  de Carvalho, M. and Martos, G. (2018). Brexit: Tracking and
  disentangling the sentiment towards leaving the EU. Submitted.
}
\examples{
## SIMULATED DATA EXAMPLES
## Example 1
N <- 500
t <- seq(.01, 5, length = N)
Y <- t^3 - 9 * t^2 + 23 * t + rnorm(N, 0, 1)
y <- mtsframe(date = t, Y)
fit <- sst(y)
plot(t, Y, col = "gray", ylab = "", 
     xlab = "time", pch = 16, ylim = c(-10, 31))
lines(t, t^3 - 9 * t^2 + 23 * t)
lines(t, fit$trendline$Y, col = "red")
## Example 2
N <- 500
t <- seq(.01, 5, length = N)
Y <- 10 * sin(3 * t) / t + rnorm(N, 0, 1)
y <- mtsframe(date = t, Y)
fit <- msst(y)
fit$selected.components
plot(t, Y, col = "gray", ylab = "",
     xlab = "time", pch = 16, ylim = c(-10, 31))
lines(t, 10 * sin(3 * t) / t, col = "black")
lines(t, fit$trendline$Y, col = "red")

## BREXIT DATA EXAMPLE
## Note: sst also can deal with several time series as input yet the
## most appropriate method for multivariate time series is msst
data(brexit)
attach(brexit)
y <- mtsframe(date, brexit[, 1:3] / 100)
fit <- sst(y)
plot(fit)

## Number of components automatically selected 
fit$selected.components 

## Chronological plot (de Carvalho and Martos, 2018; Fig. 1)
plot(fit, options = list(type = "trendlines"), xlab = "time",
     col = c("blue", "red", "black"), lwd = 2, lty = c(1, 2, 3))

## Plot elementary reconstructed components
## (de Carvalho and Martos, 2018; Fig. 3)
plot(fit, options = list(type = "components", serie = 1, ncomp = 1:2))
plot(fit, options = list(type = "components", serie = 2, ncomp = 1:2))
plot(fit, options = list(type = "components", serie = 3, ncomp = 1  ))

## Plot cumulative periodograms (with 95% confidence bands)
par(mfrow=c(1,3))
plot(fit, options = list(type = "cpgrams",
           series.names = c('Leave','Stay','Undecided')))

## Scree-plot (with 95% confidence bands)
par(mfrow = c(1, 3))
plot(fit, options = list(type = "screeplots", ncomp = 1:10,
          series.names = c('Leave','Stay','Undecided')),
            type = "b", pch = 20, lwd = 2)
}
\keyword{datasets}
\author{Gabriel Martos and Miguel de Carvalho}