\encoding{UTF-8}
\name{bmssa}
\alias{bmssa}
\alias{bmssa.default}
\title{Multivariate Singular Spectrum Business Cycle Indicator}
\description{
  Computes a business cycle indicator using multivariate singular
  spectrum analysis.
}
\usage{
bmssa(y, l = 32)
}
\arguments{
  \item{y}{multivariate time series of economic activity data from which
    the cycle is to be extracted; the first column is reserved to Gross
    Domestic Product (GDP).}
  \item{l}{window length; by default, \code{l = 32}.}
}
\value{
  \item{cycle}{time series with the business cycle indicator.}
  \item{sfisher}{vector with indices of elementary reconstructed
    components selected with Fisher \eqn{g} statistic; see details.}
  \item{erc}{time series with elementary reconstructed components
    resulting from targeted grouping based on a Fisher \eqn{g}
    statistic.}
  \item{l}{window length.}
}
\details{
  The business cycle indicator produced using this routine is based on
  methods proposed in de Carvalho and Rua (2017). A quick summary of the
  method is as follows. Multivariate singular spectrum analysis is used
  to decompose a multivariate time series (\code{y}) into principal
  components, and a Fisher \eqn{g} statistic automatically selects
  elementary reconstructed components (\code{erc}) within business cycle
  frequencies. The indicator results from adding elementary
  reconstructed components within business cycle frequencies. The
  \code{plot} method depicts the resulting business cycle indicator,
  and the \code{print} method reports the business cycle indicator
  along with the components selected by the Fisher \eqn{g} statistic.
}
\references{
  de Carvalho, M., Rodrigues, P., and Rua, A. (2012). Tracking the US
  business cycle with a singular spectrum analysis. \emph{Economics
  Letters}, \bold{114}, 32--35.

  de Carvalho, M. and Rua, A. (2017). Real-time nowcasting the US output
  gap: Singular spectrum analysis at work. \emph{International Journal
  of Forecasting}, \bold{33}, 185--198.
}
\seealso{
  See \code{\link{combplot}} for a chart of the selected elementary
  reconstructed components from which the business cycle indicator
  results.  See \code{\link{bssa}} for a univariate version of the
  method.
}
\examples{
## Tracking the US Business Cycle (de Carvalho et al, 2017; Fig. 6) 
data(GDPIP)
fit <- bmssa(log(GDPIP))
plot(fit)
print(fit)
}
\author{Miguel de Carvalho.}