\encoding{UTF-8}
\name{mtsframe}
\alias{mtsframe}
\alias{mtsframe.default}
\title{Multivariate Time Series Frame Objects}
\description{
  The function \code{mtsframe} is used to create mutivariate time series
  objects to be used in combination with the functions in the package
  ASSA.
}
\usage{
mtsframe(dates, Y)
}
\arguments{
  \item{dates}{dates at which observations took place (in case of using
    an index, the format of this field will be taken as numeric.)}
  \item{Y}{matrix with time-series in columns and observations in rows.}
}
\examples{
data(brexit)
attach(brexit)
y <- mtsframe(date, Y = brexit[, 1:3])
print(y)
# Ploting the time series under study (blue = Leave, red = Stay, black = Undecided)
plot(y, col = c('blue', 'red', 'black'), time.format = '\%m-\%y')
# When 'date' is a time index, the time.format specification is omitted.
}
\keyword{multivariate time series formatting.}