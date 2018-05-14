\encoding{UTF-8}
\name{combplot}
\alias{combplot}
\title{Comb-plot}
\description{
  Produces a comb-plot for visualizing what principal components are
  used for producing the (multivariate) singular spectrum business cycle
  indicator.
}
\usage{
combplot(fit)
}
\arguments{
  \item{fit}{a \code{bssa} or a \code{bmssa} object.}
}
\details{
  \code{combplot} yields a comb-plot indentifying the indices 
  of the components selected according to the Fisher \eqn{g} statistic,
  along with the corresponding principal components; see de Carvalho and
  Rua (2017, p. 190) for a definition.
}
\references{
  de Carvalho, M. and Rua, A. (2017). Real-time nowcasting the US output
  gap: Singular spectrum analysis at work. \emph{International Journal
  of Forecasting}, \bold{33}, 185--198.
}
\seealso{
  \code{\link{bssa}}.
}
\keyword{
  comb-plot
}
\examples{
## Tracking the US Business Cycle (de Carvalho and Rua, 2017; Fig. 5)
data(GDPIP)
fit <- bssa(log(GDPIP[, 1]))
combplot(fit)
}
\author{Miguel de Carvalho.}