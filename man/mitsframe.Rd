\encoding{UTF-8}
\name{mitsframe}
\alias{mitsframe}
\alias{mitsframe.default}
\title{Multivariate Interval Valued Time Series Frame Objects}
\description{
  The function \code{mitsframe} is used to create interval-valued multivariate time series objects.
}
\usage{
mitsframe(dates, A, B)
}
\arguments{
  \item{dates}{a vector of dates at which observations took place.} 
  \item{A}{matrix with time series in columns (dimensions) and observations in rows corresponding to the lowest values. Data must be sorted in ascendant way (first row in 'A' corresponds to the oldest values of the series and last row in 'A' corresponds to the newest values).}
  \item{B}{matrix with time series in columns (dimensions) and observations in rows corresponding to the lowest values. Data must be sorted in ascendant way (first row in 'A' corresponds to the oldest values of the series and last row in 'A' corresponds to the newest values).}
}
\references{
de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
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

plot(y) # standard plot.
}
\author{Gabriel Martos and Miguel de Carvalho.}
