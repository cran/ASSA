\encoding{UTF-8}
\name{mtsframe}
\alias{mtsframe}
\alias{mtsframe.default}
\title{Multivariate Time Series Frame Objects}
\description{
  The function \code{mtsframe} create a mutivariate time series object to be used in combination with the functions in the package ASSA.
}
\usage{
mtsframe(dates, Y)
}
\arguments{
  \item{dates}{dates at which observations took place.}
  \item{Y}{matrix with different time-series in columns (dimensions) and observations in rows. Values must be be sorted in ascendant way (first row in 'Y' corresponds to the oldest values of the series and last row in 'Y' corresponds to the newest values).}
}
\examples{
data(brexit); attach(brexit)
head(brexit, 3)
y <- mtsframe(date, Y = brexit[, 1:3])
print(y) # A 'list' with 4 elements: dates, series data matrix, and series length.

head(y$Y, 3)
y$n
y$D

plot(y) # standard plot.

# Customized plot (time.format is an additional feature to use when y$date is in 'date' format)
plot(y, time.format = '\%Y' , col = c('blue','red','black'), lty = 2, type = 'p',
     pch = 20, main = 'Brexit data', xlab = 'Year', ylab ='Trendline estimations')
}
\keyword{multivariate time series formatting.}