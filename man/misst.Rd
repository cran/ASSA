\encoding{UTF-8}
\name{misst}
\alias{misst}
\alias{misst.default}
\title{Multivariate Interval Singular Spectrum Trendlines}
\description{
  Computes a trendline for multivariate interval data using singular
  spectrum analysis.
}
\usage{
misst(y, l= 'automatic' , m = 'automatic', vertical = TRUE)
}
\arguments{
  \item{y}{object of class \code{mitsframe} (multivariate interval time series data).} 
  \item{l}{window length; the string \code{"automatic"} sets
    the default option \cr \code{l = ceiling(y$n + 1) / (y$D+1)} for
    \code{vertical} and \code{ceiling(D * (y$n + 1) / (y$D + 1))}.}
  \item{m}{number of leading eigentriples. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default by using the string \code{"automatic"}.}
  \item{vertical}{logical; if \code{TRUE} the trajectory matrices are
    stacked vertically, otherwise the bind is horizontal.}
}
\details{
  Multivariate singular spectrum analysis is used to decompose interval time
  series data (\code{y}) into principal components, and a cumulative
  periodogram-based criterion automatically learns about what elementary
  reconstructed components (\code{erc}) contribute to the signal; see de
  Carvalho and Martos (2018) for details. The trendline results from
  adding elementary reconstructed components selected by the cumulative
  periodogram of the residuals. The \code{plot} method depicts the
  trendlines, and the \code{print} method reports the trendlines along
  with the components selected by the cumulative periodogram-based
  criterion.
}	
\value{
\item{trendline}{mitsframe object with interval trendline estimation from targeted grouping based on a cumulative periodogram criterion (or according to the number of components specified in vector \code{m}).}
\item{l}{window length.}
\item{m}{vector with number of components selected on each dimension.}
\item{vertical}{flag indicating if the trajectory matrices where stacked vertically.}
\item{residuals}{mitsframe object with the interval residuals from targeted grouping based on a cumulative periodogram criterion (or according to the number of components specified in vector \code{m}).}
\item{svd}{the Singular Value Decomposition of the trajectory matrix.}
\item{erc}{list with elementary reconstructed components.}
\item{observations}{mitsframe object with the raw data \code{y}.}
}
\references{
de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
}
\seealso{
  See \code{\link{msst}} for a similar routine yielding trendlines for
  standard multivariate time series of data.
}
\examples{
muX.a = function(t){ 8 + t + sin(pi*t) }    ;    muX.b = function(t){ muX.a(t) + 2 } 
muY.a = function(t){sqrt(t) + cos(pi*t/2) }    ;  muY.b = function(t){ 2*muY.a(t) + 2 }
N = 100; t=seq(0.1,2*pi,length = N);
set.seed(1)
e.x = rnorm(100); e.y = rnorm(100);
a.X = muX.a(t) + e.x; b.X = a.X + 2
a.Y = muY.a(t) + e.y        ; b.Y = 2*a.Y + 2

A <- cbind(a.X, a.Y); B <- cbind(b.X, b.Y)
y <- mitsframe(dates=t, A=A, B = B)

fit <- misst(y)
fit$l; 
fit$m; 
fit$vertical

# Estimated trendlines:
head(fit$trendlines$A,5) 
head(fit$trendlines$B,5) 

## Estimated interval trendlines
plot(fit)

## Scree-plot
plot(fit, options = list(type = "screeplots"))

## Per
plot(fit, options = list(type = "cpgrams"))

## ERC
plot(fit, options=list(type='components',ncomp=1:3))



##################################
### Forecasting with misst     ###
##################################
pred = predict(fit, p = 5)
pred$forecasts # Forecast organized in an array.
# End
}
\keyword{datasets}