\encoding{UTF-8}
\name{isst}
\alias{isst}
\alias{isst.default}
\title{Interval Singular Spectrum Trendlines}
\description{
  Computes the trendline estimation for interval time series data using singular
  spectrum analysis.
}
\usage{
isst(y, l= 'automatic' , m = 'automatic')
}
\arguments{
\item{y}{\code{itsmframe} data corresponding to univariate interval time series data.} 
\item{l}{{window length; the string \code{'automatic'} sets the
      default option \cr \code{l = ceiling(y$n + 1) / 2}.}}
\item{m}{number of leading eigentriples; the string \code{"automatic"}
    yields an automatic criterion for choosing m based on the cumulative periodogram of the
    residuals; see details.}
}
\value{
  \item{trendline}{itsframe object with interval trendline estimation from
    targeted grouping based on a cumulative periodogram criterion (or
    according to the number of components specified in \code{m}).}
 \item{l}{window length.}
 \item{m}{number of leading eigentriples. An automatic
    criterion based on the cumulative periodogram of the residuals is
    provided by default by using the string \code{"automatic"}.}
 \item{residuals}{itsframe object with the residuals from targeted
    grouping based on a cumulative periodogram criterion (or according
    to the number of components specified in \code{m}).}
  \item{svd}{Singular value decomposition corresponding to the trajectory matrix.}
  \item{erc}{elementary reconstructed components.}
 \item{observations}{itsframe object with the raw data.}
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
  See \code{\link{misst}} for a version of the routine for multivariate
  interval value time series.
} 
\references{
de Carvalho, M. and Martos, G. (2020). Modeling Interval Trendlines: Symbolic Singular Spectrum Analysis for Interval Time Series. Submitted (available on arXiv).
}
\examples{
# Merval data example:
data(merval)
id.data <-  rev(which(merval[,1]>'2015-12-31' & merval[,1]<'2020-10-01') ) 
y <-  itsframe(date=merval[id.data,1], a=merval[id.data,2], b=merval[id.data,3]); 

isst_output <- isst(y ,l = 'automatic', m = 'automatic')
print(isst_output)

# Estimated trendlines 
plot(isst_output)

## Scree-plot
plot(isst_output, options = list(type = "screeplot", ncomp = 1:10),
     type = "b", pch = 20, lwd = 2)

# Elementary reconstructed components
plot(isst_output, options=list(type='components',ncomp=1:3),
     xlab='Time')

# cpgram's ('a=low' and 'b=high')
plot(isst_output, options = list(type='cpgram'))
# Setting m = 'automatic' (default option) to obtain cpgrams inside the bandwiths.

##################################
###  Forecasting with isst     ###
##################################
pred <- predict(isst_output, p = 5)
head(pred$forecast,3) # Forecasted interval data.

attributes(pred)
pred$coefficients[1:5] # linear recurrence parameters.
# End
}
\keyword{datasets}