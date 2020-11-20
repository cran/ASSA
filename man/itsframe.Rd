\encoding{UTF-8}
\name{itsframe}
\alias{itsframe}
\alias{itsframe.default}
\title{Interval Time Series Frame Objects}
\description{
  The function \code{itsframe} creates a univariate interval time series
  object to be used in combination with the functions in the package
  ASSA.
}
\usage{
itsframe(dates, a, b)
}
\arguments{
  \item{dates}{dates at which observations took place.}
  \item{a}{vector with lower interval time-series values sorted in ascendant way (first element in 'a' corresponds to the oldest value of the interval series and last element in 'a' corresponds to the newest value).}
  \item{b}{vector with upper interval time-series values sorted in ascendant way (first element in 'b' corresponds to the oldest value of the interval series and last element in 'b' corresponds to the newest value).}
}
\references{
de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
}
\examples{
data(merval)
id.data <-  rev(which(merval[,1]>'2015-12-31' & merval[,1]<'2020-10-01') ) 
y <-  itsframe(date=merval[id.data,1], a=merval[id.data,2], b=merval[id.data,3]); 

plot(y, main = 'MERVAL')
}
\keyword{univariate interval time series formatting.}
