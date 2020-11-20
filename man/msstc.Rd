\encoding{UTF-8}
\name{msstc}
\alias{msstc}
\alias{msstc.default}
\title{Multivariate Singular Spectrum Trendlines for Compositional Data}
\description{
  Computes trendlines on the unit simplex for multivariate time series
  data using multivariate singular spectrum analysis.
}
\usage{
msstc(y, l = 'automatic', m = 'automatic', vertical = TRUE)
}
\arguments{
  \item{y}{mtsframe object containing data.}
  \item{l}{window length; the string \code{"automatic"} sets
    the default option \cr \code{l = ceiling((y$n + 1) / y$D)} for
    \code{vertical} and \code{ceiling(y$D * (y$n + 1) / (y$D + 1))} in the 
     case of horizontal binding.}
 \item{m}{vector with the number of leading eigentriples on each dimension. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default by using the string \code{"automatic"}, see details.}
  \item{vertical}{logical; if \code{TRUE} the trajectory matrices are
    stacked vertically, otherwise the bind is horizontal.}
}
\details{
  The trendline produced using this routine is based on the methods
  proposed in de Carvalho and Martos (2018). A quick summary of the
  method is as follows. Multivariate singular spectrum analysis is
  used to decompose time series data (\code{y}) into principal
  components, and a cumulative periodogram-based criterion
  automatically learns about what elementary reconstructed components
  (\code{erc}) contribute to the signal; see de Carvalho and Martos
  (2018) for details. The trendline results from adding elementary
  reconstructed components selected by the cumulative periodogram, and
  after projecting into the unit simplex. The \code{plot} method
  depicts the trendlines, and the \code{print} method reports the
  trendlines along with the components selected by the cumulative
  periodogram-based criterion.
}
\value{
\item{trendline}{mtsframe object with trendline estimation from targeted grouping based on a cumulative periodogram criterion (or according to the number of components specified in vector \code{m}).}
\item{l}{window length.}
\item{m}{vector with number of components selected on each dimension.}
\item{vertical}{flag indicating if the trajectory matrices where stacked vertically.}
\item{residuals}{mtsframe object with the residuals from targeted grouping based on a cumulative periodogram criterion (or according to the number of components specified in vector \code{m}).}
\item{svd}{the Singular Value Decomposition of the trajectory matrix.}
\item{erc}{list with elementary reconstructed components.}
\item{observations}{mtsframe object with the observations \code{y}.}
}
\references{
  de Carvalho, M. and Martos, G. (2020). Brexit: Tracking and
  disentangling the sentiment towards leaving the EU. \emph{International Journal
  of Forecasting}, \bold{36}, 1128--1137.
}
\seealso{
  See \code{\link{msst}} for a similar routine yielding trendlines for
  multivariate time series, but which does not project the pointwise
  estimates to the unit simplex.
}
\examples{
## BREXIT DATA EXAMPLE
## (de Carvalho and Martos, 2018; Fig. 1)
data(brexit)
attach(brexit)
y <- mtsframe(date, brexit[, 1:3] / 100)
fit <- msstc(y)
# Estimations on the simplex
rowSums(fit$trendlines$Y)

# Forecast also in the simplex
rowSums(predict(fit, p = 5)$forecast)


## Window length and components automatically selected
fit$l; fit$m

## Plot trendlines (de Carvalho and Martos, 2018; Fig. 1)
plot(fit, options = list(type = "trendlines"), xlab="time",
     col=c("blue", "red", "black"), lwd = 2, lty = c(1, 2, 3))


## Plot cumulative periodograms (with 95% confidence bands)
par(mfrow = c(1, 3))
plot(fit, options = list(type = "cpgrams") )

## Scree-plot (with 95% confidence bands)
par(mfrow = c(1, 1))
plot(fit, options = list(type = "screeplot", ncomp.scree = 1:10),
     type = "b", pch = 20, lwd = 2, main='Scree plot')

## Plot elementary reconstructed components 
## (de Carvalho and Martos, 2020; Fig. 5)
plot(fit, options = list(type = "components", ncomp = 1:2))
}
\keyword{datasets}
\author{Gabriel Martos and Miguel de Carvalho}
