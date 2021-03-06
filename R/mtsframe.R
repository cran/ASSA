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

plot.mtsframe <- function(x, time.format="%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                         ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL,cex.lab=NULL, cex.axis=NULL,cex.main=NULL, ...) {
    matplot(x$dates, x$Y,
             main = if(is.null(main)){''} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){'Time'} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c(min(x$Y), max(x$Y)) } else {ylim},
             xlim = if(is.null(xlim)){c(min(x$dates),max(x$dates))} else {xlim},  
             cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main},
             xaxt="n")
     if(inherits(x$dates, "Date")==T){  timelabels<-format(x$dates,time.format) ; axis(1,at=x$dates,labels=timelabels,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$dates,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
}
