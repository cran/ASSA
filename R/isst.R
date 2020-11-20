isst <-function(y, l = 'automatic', m = 'automatic')
  UseMethod("isst")

isst.default <- function(y, l = 'automatic' , m = 'automatic') {
  
n <- y$n; D <- 1
  
## Run a basic input validation
if (class(y) %notin% c('itsframe') ){
    stop('y must be a univariate interval time series object.')}
  if(is.numeric(m) & length(m) != 1) { stop('Stop: Parameter m in the model must be a postive integer.') }
  if(!is.numeric(l)){  l <- ceiling( (n + 1) / 2) }
  k <- n - l + 1
  
## Step 1: Trajectory matrix *array*
  Y <- array(NA, c(l, k, 2) ) # l rows, k columns, 2 dimensions
  Y[ , , 1] <- trajectory(y$a, l, k)
  Y[ , , 2] <- trajectory(y$b, l, k) 
  ## Step 2:svd
S = Var.Est(Y,Y)$outer.product #*array input*
SVD <- svd(S)
rk <- qr(SVD$u)$rank # rank

## Step 3 and 4:
y.tilde = residuals = matrix(0, ncol = 2, nrow = n)
if(is.numeric(m)) {
  if(m > rk) { m = rk
    warning( "Number of components (see 'm') required to construct trendlines 
              was reduced to be compatible with the rank of the trajectory matrix." ) } 
      Y_i = YYI_i(Y, SVD, i=m, l=l, k=k) ; 
      y.tilde[,1]     <- dbar(Y_i$Y_iA, l = l, k = k);  
      y.tilde[,2]     <- dbar( Y_i$Y_iB, l = l, k = k)
      residuals[,1]   <- pmin(y$a-y.tilde[,1],y$b-y.tilde[,2]);  
      residuals[,2]   <- pmax(y$a-y.tilde[,1],y$b-y.tilde[,2]);
          } else { m = 0;  stop.flag = 1
          while(stop.flag != 0 & m < rk) { 
            m = m + 1
            Y_i = YYI_i(Y, SVD, i=m, l=l, k=k) ; 
            y.tilde[,1]     <- dbar(Y_i$Y_iA, l = l, k = k);  
            y.tilde[,2]     <- dbar(Y_i$Y_iB, l = l, k = k);
            residuals[,1]   <- pmin(y$a-y.tilde[,1],y$b-y.tilde[,2]);  
            residuals[,2]   <- pmax(y$a-y.tilde[,1],y$b-y.tilde[,2]);
            stop.flag <- ecip(residuals);  
          } } # end 'else'

  y.tilde   <- itsframe(dates = y$date, a = pmin(y.tilde[ , 1],y.tilde[ , 2]), b = pmax(y.tilde[, 1], y.tilde[, 2]))
  residuals <- itsframe(dates = y$date, a = residuals[ , 1], b = residuals[, 2])

#### ERC:
erc <- array(NA, c(n, m, 2) ) # n rows, m columns, 2 dimensions
for(i in 1:m){      
Y_i = YYI_i.erc(Y, SVD, i=i, l=l, k=k) ; 
ll = dbar(Y_i$Y_iA, l = l, k = k); uu = dbar( Y_i$Y_iB, l = l, k = k)
erc[,i,1] <- pmin(ll,uu); erc[,i,2] <- pmax(ll,uu)
} 
    
outputs <- list(trendline = y.tilde, 
                  l = l, 
                  m = m, 
                  residuals = residuals,
                  svd = SVD,
                  erc = erc,
                  observations = y,
                  call = match.call()) ; 
  class(outputs) <- "isst"
  return(outputs)
  
}

print.isst <- function(x, ...) {
  cat("\n Singular Spectrum Trendlines for Interval Data:\n ========================================= \n")
  print(x$call)
  cat("\n Interval Trendlines:\n"); cat("\n")
  print(x$trendline)
  cat("\n Elementary components included in the estimation \n"); cat("\n")
  print(x$m)
}

plot.isst <- function(x, time.format = "%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                      ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL, cex.lab = NULL, cex.axis = NULL, cex.main = NULL,  
                      options = list(type = 'trendline', 
                                     ncomp = NULL), ...) {
  if(options$type %!in% c('trendline', 'screeplot', 'components', 'cpgram')) {
    stop('options type must be one of the strings: "trendline", "components", "cpgram", or "screeplot".')}
  
  if(options$type == 'trendline') {
    plot(x$trendline$date, x$trendline$a,
         main = if(is.null(main)){''} else {main},
         col  = if(is.null(col)){'black'} else {col},
         lty  = if(is.null(lty)){1} else {lty},
         pch  = if(is.null(pch)){1} else {pch},
         lwd  = if(is.null(lwd)){1} else {lwd},
         type = if(is.null(type)){'l'} else {type},
         xlab = if(is.null(xlab)){'Time'} else {xlab},
         ylab = if(is.null(ylab)){'Interval Singular Spectrum Trendline'} else {ylab}, 
         ylim = if(is.null(ylim)){c(min(x$trendline$a,x$trendline$b),max(x$trendline$a,x$trendline$b))} else {ylim}, 
         xlim = if(is.null(xlim)){c(min(x$trendline$dates),max(x$trendline$dates))} else {xlim}, 
         cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
         cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
         cex.main = if(is.null(cex.main)){1} else {cex.main},
         xaxt = "n")
    if(inherits(x$trendline$dates, "Date")==T){  timelabels<-format(x$trendline$dates,time.format) ;    axis(1,at=x$trendline$dates,labels=timelabels,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$trendline$dates,cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
    
    lines(x$trendline$date, x$trendline$b, 
          col  = if(is.null(col)){'black'} else {col},
          lty  = if(is.null(lty)){1} else {lty},
          pch  = if(is.null(pch)){1} else {pch},
          type = if(is.null(type)){'l'} else {type})
    
    polygon(c(rev(x$trendline$date), x$trendline$date), c(rev(x$trendline$a), x$trendline$b), 
            border = NA, col = if(is.null(col)){'lightgray'} else {col})   
  }
  
  
  
  if(options$type == 'screeplot') {
    if(is.null(options$ncomp)){options$ncomp = c(1:3)}
    plot(log(x$svd$d[options$ncomp]), 
         main = if(is.null(main)){'Scree--plot'} else {main},
         col  = if(is.null(col)){'black'} else {col},
         lty  = if(is.null(lty)){1} else {lty},
         pch  = if(is.null(pch)){1} else {pch},
         lwd  = if(is.null(lwd)){1} else {lwd},
         type = if(is.null(type)){'b'} else {type},
         cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
         cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
         cex.main = if(is.null(cex.main)){1} else {cex.main},
         xlab = if(is.null(xlab)){'Index'} else {xlab},
         ylab = if(is.null(ylab)){'Eigenvalues (log-scale)'} else {ylab})
  } 
  
  if(options$type == 'components') {
    if(is.null(options$ncomp)){ options$ncomp = c(1:round(dim(x$erc)[2])) }
    if( max(options$ncomp) > dim(x$erc)[2] ) {
      print( paste('Please choose the number of ERC in the range from 1 :', dim(x$erc)[2]) ) 
    } else {
      coldim <- options$ncomp;     
      par(mfrow = c(1, length(coldim)))
      for(j in coldim) {
        plot(x$trendline$dates, x$erc[, j, 1],
             main = if(is.null(main)){paste("ERC -", j)} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){''} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c(min(c(x$erc[, j, 1],x$erc[, j, 2])),max(c(x$erc[, j, 1],x$erc[, j, 2]))) } else {ylim},
             cex.lab = if(is.null(cex.lab)){0.8} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){0.8} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main} )
        lines(x$trendline$dates, x$erc[, j, 2],  
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type})
        polygon(c(rev(x$trendline$dates), x$trendline$dates),
                c(rev(x$erc[, j, 1]), x$erc[, j, 2]), 
                border = NA, col = if(is.null(col)){'lightgray'} else {col})
      }
    }
 }
    
    if(options$type == "cpgram") {
      ecip(residuos = cbind(x$residuals$a,x$residuals$b), plot.flag = TRUE)
    }    
  }
