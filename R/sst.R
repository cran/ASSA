sst <- function(y, l = 'automatic', m = 'automatic')
  UseMethod("sst")

sst.default <- function(y, l = 'automatic', m = 'automatic') {  
  n <- y$n
  if (class(y) != "tsframe") { stop('Input data must be an tsframe object.') }
  if(!is.numeric(l)){  l <- ceiling( (n + 1) / 2) }
  k <- n - l + 1
  if(!is.numeric(m) & m !='automatic'){ stop('Plase choose the automatic criterion or an integer value for m.')} 
  if(is.numeric(m) & (length(m) != 1 | m<0  ) ) { stop('m must be positive integer.') }
  ## Step 1:
  Y <- trajectory(y$y, l, k)

  ## Step 2:
  SVD <- svd(Y%*%t(Y))

  ## Step 3 and 4:
  erc = list(); # *will contain one element (a matrix with the ERC by rows)*    
  if(is.numeric(m)) {
     erc[[1]] <- matrix(NA, nrow = n,  ncol = m) 
     for(i in 1:m){ y_i <- Y_i(Y = Y, SVD = SVD, i = i); erc[[1]][,i] <- dbar(y_i, l = l, k = k) }
    y.tilde   <- rowSums(erc[[1]])
    residuals <- y.tilde - y$y # residuals.
           } else {  
              stop.flag <- 0; mm <- 1; erc[[1]] <- matrix(NA, nrow = n,  ncol = 0)
              while(stop.flag == 0) {
              y_i <- Y_i(Y = Y, SVD = SVD, i = mm); erc[[1]] <-cbind(erc[[1]], dbar(y_i, l = l, k = k))  
              y.tilde   <- rowSums(erc[[1]])
              residuals <- y.tilde - y$y # residuals.
              if(cpgram2(residuals) == 1) { mm <- mm + 1 } else { stop.flag <- 1 }
                              } # end 'while'
              m = mm } # end 'else'
  
  y.tilde <- tsframe(dates = y$date, y = y.tilde)
  residuals <- tsframe(dates = y$date, y = residuals)
  outputs <- list(trendline = y.tilde, 
                  l = l, 
                  m = m, 
                  residuals = residuals,
		  svd = SVD,
 		  erc = erc,
                  observations = y,
                  call = match.call()) ; 
  class(outputs) <- "sst"
  return(outputs)
} 

print.sst <- function(x, ...) {
  cat("\n Singular Spectrum Trendline:\n ============================ \n")
  print(x$call)
  cat("\n Trendlines:\n"); cat("\n")
  print(x$trendline)
  cat("\n Elementary components included in the estimation \n"); cat("\n")
  print(x$m)
}

plot.sst <- function(x, time.format="%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                        ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL, cex.lab = NULL, cex.axis = NULL, cex.main = NULL,
                        options = list(type = "trendline", ncomp = NULL), ...) {
# Input check:
  if(options$type %!in% c("trendline","components","screeplot","cpgram")){ 
    stop("options type must be one of the strings: 'trendline', 'components', 'cpgram', or 'screeplot'.") }

  if(options$type=="trendline") {    
   plot(x$trendline$dates, x$trendline$y, 
              main = if(is.null(main)){''} else {main},
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              lwd  = if(is.null(lwd)){1} else {lwd},
              type = if(is.null(type)){'b'} else {type},
              xlab = if(is.null(xlab)){'Time'} else {xlab},
              ylab = if(is.null(ylab)){'Singular Spectrum Trendline'} else {ylab}, 
              ylim = if(is.null(ylim)){c(min(x$trendline$y),max(x$trendline$y))} else {ylim}, 
              xlim = if(is.null(xlim)){c(min(x$trendline$dates),max(x$trendline$dates))} else {xlim}, 
           cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
          cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
          cex.main = if(is.null(cex.main)){1} else {cex.main},
              xaxt = "n")
  if(inherits(x$trendline$dates, "Date")==T){  timelabels<-format(x$trendline$dates,time.format) ;     axis(1,at=x$trendline$dates,labels=timelabels)} else {  axis(1,at=x$trendline$dates)    }
  }
  
  if(options$type == 'screeplot') {
      if(is.null(options$ncomp)){options$ncomp = (1:3)}
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
        if(is.null(options$ncomp)){options$ncomp = c(1:round(dim(x$erc[[1]])[2]))}
        if( max(options$ncomp) > dim(x$erc[[1]])[2] ) {
            print( paste('Please choose the number of ERC in the range from 1 :', dim(x$erc[[1]])[2]) ) 
          } else {
           coldim <- options$ncomp
           par(mfrow = c(1, length(coldim)))
                for(j in coldim) {
                   plot(x$erc[[1]][,j],
              main = if(is.null(main)){paste("ERC -", j)} else {main},
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              lwd  = if(is.null(lwd)){1} else {lwd},
              type = if(is.null(type)){'l'} else {type},
              xlab = if(is.null(xlab)){''} else {xlab},
              ylab = if(is.null(ylab)){''} else {ylab},
           cex.lab = if(is.null(cex.lab)){0.8} else {cex.lab},
          cex.axis = if(is.null(cex.axis)){0.8} else {cex.axis},
          cex.main = if(is.null(cex.main)){1} else {cex.main} )
           }
      }
  }
 
   
  if(options$type == "cpgram") { cpgram(x$residuals$y,main = "") }
}
