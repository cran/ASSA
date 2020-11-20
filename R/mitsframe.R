mitsframe <- function(dates, A,B)
  UseMethod("mitsframe")

mitsframe.default <- function(dates, A,B) {
  A = as.matrix(A) # force to be a matrix in case of univariate interval time series.
  B = as.matrix(B) # force to be to a matrix in case of univariate interval time series.
  n <- dim(A)[1]; D <- dim(A)[2]
  
  if(inherits(dates, "Date")==T){dates = as.Date(dates)} else {dates = dates}
  ## Organize and return outputs
  outputs <- list(dates = dates, A = A, B=B,  n=n, D=D, call = match.call())
  class(outputs) <- "mitsframe"
  return(outputs)
}

plot.mitsframe <- function(x,time.format="%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                           tick = TRUE, ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL,cex.lab=NULL, cex.axis=NULL,cex.main=NULL, ...) {
      matplot(x$dates, x$A,
             main = if(is.null(main)){''} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){'Time'} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c(min(x$A,x$B), max(x$A,x$B)) } else {ylim},
             xlim = if(is.null(xlim)){c(min(x$dates),max(x$dates))} else {xlim},  
             cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main},
             xaxt="n")
     if(inherits(x$dates, "Date")==T){  timelabels<-format(x$dates,time.format) ; 
     axis(1,at=x$dates,tick =tick, labels=timelabels,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$dates,tick =tick, cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
  
  for(i in 1:x$D){
  lines(x$dates, x$B[,i],
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              type = if(is.null(type)){'l'} else {type},
              lwd = if(is.null(lwd)){1} else {lwd})
  polygon(c(rev(x$dates), x$dates), c(rev(x$A[,i]), x$B[,i]), 
          border = NA, col = if(is.null(col)){'lightgray'} else {col})
  }
}
