\encoding{UTF-8}
\name{merval}
\alias{merval}
\title{MERVAL interval data}
\description{
Raw interval data series corresponding to weekly minimum and maximum values of MERVAL index ranging from January 1st 2016 to September 30th 2020.
}
\usage{merval}
\format{
  A dataframe with 353 observations. The object
  is of class \code{list}.}
\references{
de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
}
\source{
Yahoo Finance. 
}	
\examples{
data(merval)
head(merval,3)
}
\keyword{datasets}