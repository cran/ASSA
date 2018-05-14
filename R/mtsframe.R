mtsframe <- function(dates, Y)
    UseMethod("mtsframe")

mtsframe.default <- function(dates, Y) {
  Y = as.matrix(Y) # coherece to a matrix in case of univariate time series. 
  n <- dim(Y)[1]; D <- dim(Y)[2] 
  if(inherits(dates, "Date")==T){dates = as.Date(dates)} else {dates = dates}
      ## Organize and return outputs
    outputs <- list(dates = dates, Y = Y, n=n, D=D, call = match.call())
    class(outputs) <- "mtsframe"
    return(outputs)
}

plot.mtsframe <- function(x, xlab = "Time", ylab = "", time.format="%m-%y", ...) {
    matplot(x$dates, x$Y, xlab = "Time",ylab = "", type = "l", 
         ylim = c(min(x$Y), max(x$Y)),xaxt="n",...)
  
  if(inherits(x$dates, "Date")==T){
    timelabels<-format(x$dates,time.format)
    axis(1,at=x$dates,labels=timelabels)} else {
      axis(1,at=x$dates)    }
}