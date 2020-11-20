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
  \item{y}{tsframe format data containing univariate time series data. 
           More appropriate method for multivariate time series is \code{\link{msst}}. }
  \item{l}{window length; the string \code{"automatic"} automatic sets
    the default option \cr \code{l = ceiling(y$n + 1) / 2}.}
  \item{m}{number of leading eigentriples; the string \code{"automatic"}
    yields an automatic criterion for choosing m based on the cumulative periodogram of the
    residuals; see details.}}
\value{
  \item{trendline}{tsframe object with trendline estimation from
    targeted grouping based on a cumulative periodogram criterion (or
    according to the number of components specified in \code{m}).}
 \item{l}{window length.}
 \item{m}{number of leading eigentriples. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default by using the string \code{"automatic"}.}
 \item{residuals}{tsframe object with the residuals from targeted
    grouping based on a cumulative periodogram criterion (or according
    to the number of components specified in \code{m}).}
  \item{svd}{Singular value decomposition corresponding to the trajectory matrix.}
  \item{erc}{elementary reconstructed components.}
 \item{observations}{tsframe object with the raw data.}
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
  a list including the strings \code{"trendline"},
  \code{"components"}, \code{"cpgram"}, and \code{"screeplot"},
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
   de Carvalho, M. and Martos, G. (2020). Brexit: Tracking and
  disentangling the sentiment towards leaving the EU. \emph{International Journal
  of Forecasting}, \bold{36}, 1128--1137.
}
\examples{
## BREXIT DATA EXAMPLE
data(brexit); attach(brexit)
l <- tsframe(date, brexit[, 1] / 100) # l = leave
fit <- sst(l); 
fit$m; fit$l # Number of ERC and parameter l in the model.
plot(fit, col = "red", lwd = 3, xlab = 'Time', ylab = 'Leave')
points(date, brexit[, 1] / 100, pch = 20)

## Scree-plot
plot(fit, options = list(type = "screeplot", ncomp = 1:10,
                         series.names = c('Leave')), type = "b", pch = 20, lwd = 2)


## Plot cumulative periodogram
par(mfrow=c(1,1), mar=c(4,2,1,1))
plot(fit, options = list(type = "cpgram", series.names = c('Leave')) )

## Elementary Reconstructed Components (ERC) plot:
plot(fit, options = list(type = "components", ncomp = 1:2))
}
\keyword{datasets}
\author{Gabriel Martos and Miguel de Carvalho}
