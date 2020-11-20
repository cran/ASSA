\encoding{UTF-8}
\name{tsframe}
\alias{tsframe}
\alias{tsframe.default}
\title{Time Series Frame Objects}
\description{
  The function \code{tsframe} creates a univariate time series
  object to be used in combination with the functions in the package
  ASSA.
}
\usage{
tsframe(dates, y)
}
\arguments{
  \item{dates}{dates at which observations took place. }
  \item{y}{vector with time-series values sorted in ascendant way (first element in 'y' corresponds to the oldest value of the series and last element in 'y' corresponds to the newest value).}
}
\examples{
data(brexit); attach(brexit)
head(brexit, 3)
y <- tsframe(date, y = brexit[, 1]) # data is  
print(y) # 'list' with 4 elements: dates, the series data, and serie length.
plot(y, col = 'blue' , lwd = 2, lty = 1)
}
\keyword{univariate time series formatting.}