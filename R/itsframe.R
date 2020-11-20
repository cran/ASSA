##  ========================================================================  ##
##  Miguel de Carvalho                                                        ##
##  Copyright (C) 2018                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

itsframe <- function(dates, a, b)
    UseMethod("itsframe")

itsframe.default <- function(dates, a, b) {
    n <- length(a)
    ## Run basic input validation
    if (length(a) != length(b))
        stop('a and b must be of the same length')
    
    if(inherits(dates, "Date")==T){dates = as.Date(dates)} else {dates = dates}
    ## Organize and return outputs
    outputs <- list(dates = dates, a = a, b = b, n=n, D=1, call = match.call())
    class(outputs) <- "itsframe"
    return(outputs)
}

plot.itsframe <- function(x, time.format="%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                          tick = TRUE, ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL,cex.lab=NULL, cex.axis=NULL,cex.main=NULL, ...) {
    plot(x$dates, x$a, 
              main = if(is.null(main)){''} else {main},
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              type = if(is.null(type)){'l'} else {type},
              xlab = if(is.null(xlab)){'Time'} else {xlab},
               lwd = if(is.null(lwd)){1} else {lwd},
              ylab = if(is.null(ylab)){''} else {ylab}, 
              ylim = if(is.null(ylim)){c(min(x$a, x$b),max(x$a, x$b))} else {ylim}, 
              xlim = if(is.null(xlim)){c(min(x$dates),max(x$dates))} else {xlim},  
           cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
          cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
          cex.main = if(is.null(cex.main)){1} else {cex.main}, 
              xaxt="n")
    if(inherits(x$dates, "Date")==T){  timelabels<-format(x$dates,time.format) ; 
    axis(1,at=x$dates, tick =tick, labels=timelabels,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$dates, tick =tick, cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
    
    lines(x$dates, x$b, 
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              type = if(is.null(type)){'l'} else {type},
              lwd = if(is.null(lwd)){1} else {lwd})
    polygon(c(rev(x$dates), x$dates), c(rev(x$a), x$b), 
            border = NA, col = if(is.null(col)){'lightgray'} else {col})   
}
