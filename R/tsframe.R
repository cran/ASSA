tsframe <- function(dates, y)
    UseMethod("tsframe")
tsframe.default <- function(dates, y) {
  y = as.matrix(y);  n <- dim(y)[1];  D <- dim(y)[2] 
 if(D != 1) { stop('Plase, give a univariate time series as an input data.') }
   
 if(inherits(dates, "Date")==T){dates = as.Date(dates)} else {dates = dates}
    outputs <- list(dates = dates, y = y, n=n, D=1, call = match.call())
    class(outputs) <- "tsframe"
    return(outputs)
}

plot.tsframe <- function(x, time.format="%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                         ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL,cex.lab=NULL, cex.axis=NULL,cex.main=NULL, ...) {
    plot(x$dates, x$y, 
             main = if(is.null(main)){''} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){'Time'} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c(min(x$y), max(x$y)) } else {ylim},
             xlim = if(is.null(xlim)){c(min(x$dates),max(x$dates))} else {xlim},  
             cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main},
             xaxt="n")
     if(inherits(x$dates, "Date")==T){  timelabels<-format(x$dates,time.format) ; axis(1,at=x$dates,labels=timelabels,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$dates,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
}
