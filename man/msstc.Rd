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
  \item{l}{{window length; the string \code{'automatic'} sets the
      default option \cr \code{l = ceiling(y$n + 1) / y$D}.}}
  \item{m}{number of leading eigentriples; the string \code{'automatic'}
    yields a vector containing the number of components in each
    dimension to be used in the trendline estimation. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default; see details.}
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
  \item{trendline}{mtsframe object with trendline estimation from
    targeted grouping based on a cumulative periodogram criterion (or
    according to the number of components specified in vector
    \code{m}).}
  \item{residuals}{mtsframe object with the residuals from targeted
    grouping based on a cumulative periodogram criterion (or according
    to the number of components specified in vector \code{m}).}
  \item{erc}{list with elementary reconstructed components.}
  \item{eigen.val}{vector with the singular values of the trajectory
    matrix.}
  \item{l}{window length.}
  \item{selected.components}{vector with number of components selected
    on each dimension.}
  \item{selection.criteria}{a vector indicating if the null hypothesis
    of white noise is rejected along the dimensions (0: not rejected,
    1: rejected).}
  \item{rank}{rank of the trajectory matrix.}
}
\references{
  de Carvalho, M. and Martos, G. (2018). Brexit: Tracking and
  disentangling the sentiment towards leaving the EU. Submitted.
}
\seealso{
  See \code{\link{msst}} for a similar routine yielding trendlines for
  multivariate time series, but which does not project the pointwise
  estimates to the unit simplex.
}
\examples{
## Brexit data and MSSA on the simplex
## (de Carvalho and Martos, 2018; Fig. 1)
data(brexit)
attach(brexit)
y <- mtsframe(date, brexit[, 1:3] / 100)
fit <- msstc(y)

## Window length and number of components automatically selected on
## each dimension:
fit$l
fit$selected.components 

## Plot trendlines (de Carvalho and Martos, 2018; Fig. 1)
plot(fit, options = list(type = "trendlines"), xlab="time",
     col=c("blue", "red", "black"), lwd = 2, lty = c(1, 2, 3))

## Plot elementary reconstructed components 
## (de Carvalho and Martos, 2018; Fig. 5)
plot(fit, options = list(type = "components", ncomp = 1:3))

## Plot cumulative periodograms (with 95% confidence bands)
par(mfrow = c(1, 3))
plot(fit, options = list(type = "cpgrams",
          series.names = c('Leave','Stay','Undecided')) )

## Scree-plot (with 95% confidence bands)
par(mfrow = c(1, 1))
plot(fit, options = list(type = "screeplots", ncomp = 1:10),
     type = "b", pch = 20, lwd = 2, main='Scree plot')

}
\keyword{datasets}
\author{Gabriel Martos and Miguel de Carvalho}