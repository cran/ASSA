msst <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE)
  UseMethod("msst")

msst.default <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE) {
  n <- y$n; D <- y$D
  ## Cf H. Hassani & R. Mahmoudvand (2013; pp. 68, eq. 25).
  if((l == 'automatic') & (vertical == TRUE))
    l <- ceiling((n + 1) / (D + 1)) 
  if((l == 'automatic') & (vertical != TRUE))
    l <- ceiling(D * (n + 1) / (D + 1)) 
  k <- n - l + 1
  
  ## Run a basic input validation
  if (class(y) != "mtsframe")
    stop('The input y must be an mtsframe object.')
  if (length(l) >1 | l%%1 != 0 | l <= 0 | l > y$n)
    stop('l is a positive  integer smaller than the observations on each series.')
  if(is.numeric(m) & length(m) != D)
    stop('Vector m must be equal to the number of time series.')
  
  ## Step 1: Embedding
  if(vertical == TRUE) {
  ## Vertical embedding
    Y <- matrix(NA, nrow = D * l, ncol = k)
    for(d in 1:D)
      Y[((d - 1) * l + 1): (d * l), ] <- trajectory(y$Y[, d], l, k)
 } else {
    ## Horizontal embedding 
    Y <- matrix(NA, nrow =  l, ncol = D * k)
    for(d in 1:D)
      Y[, ((d - 1) * k + 1) : (d * k) ] <- trajectory(y$Y[, d], l, k)
  }
  rk <- qr(Y)$rank # rank number (to make a consistency check later)
  
  ## Step 2: SVD
  SVD <- svd(Y%*%t(Y))
  
  ## Step 3 and 4:
  erc = list(); Y.tilde = Residuals = matrix(NA,nrow = n, ncol = D)
  if(is.numeric(m) == T) {
    if(sum(m > rk) > 0) { replace(m, m > rk, rk) # Rank consistency check first.
      warning("some entries in m automatically reduced--rank deficient trajectory matrix--." ) }    
    for(d in 1:D) {
      erc[[d]] <- matrix(NA, nrow = n,  ncol = m[d]) 
      for(i in 1:m[d]){ y_i <- Y_i(Y = Y, SVD = SVD, i = i); 
      if(vertical == TRUE) {erc[[d]][,i] <- dbar(y_i[((d - 1) * l + 1): (d * l), ], l = l, k = k) 
      } else {erc[[d]][,i] <- dbar(y_i[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k)}   }
      Y.tilde[,d]   <- rowSums(erc[[d]])
      Residuals[,d] <- Y.tilde[,d] - y$Y[,d] # Residuals matrix
    }
    # end if:numeric.
  } else { 
  m = c() # *store here the total number of ERC on each dimension*
    for(d in 1:D) {
      stop.flag <- 0;  mm <- 1; erc[[d]] <- matrix(NA, nrow = n,  ncol = 0) 
      while(stop.flag==0){ y_i <- Y_i(Y = Y, SVD = SVD, i = mm); 
      if(vertical == TRUE) {erc[[d]] <- cbind(erc[[d]], dbar(y_i[((d - 1) * l + 1): (d * l), ], l = l, k = k)) 
      } else {erc[[d]] <- cbind(erc[[d]],dbar(y_i[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k))}   
      Y.tilde[,d]   <- rowSums(erc[[d]])
      Residuals[,d] <- Y.tilde[,d] - y$Y[,d] # Residuals matrix
      if(cpgram2(Residuals[,d]) == 1) { mm <- mm + 1 } else { stop.flag <- 1 }
      } # end 'while'
      m[d] = mm } # end 'for:d'
  } #end steps 3 & 4.
  
  # Grouping the results into mtsframe objects to deliver in the output
  
  Y.tilde   <- mtsframe(dates = y$date, Y = Y.tilde)
  Residuals <- mtsframe(dates = y$date, Y = Residuals)
  
  outputs <- list(trendlines = Y.tilde, 
                  l = l, 
                  m = m,
                  vertical = vertical, 
                  residuals = Residuals,
                  svd = SVD,
                  erc = erc,
                  observations = y,
                  call = match.call()) ;   
  class(outputs) <- "msst"
  return(outputs)
}

print.msst <- function(x, ...) {
  cat("\n Multivariate Singular Spectrum Trendlines:\n ========================================= \n")
  print(x$call)
  cat("\n Trendlines:\n"); cat("\n")
  print(x$trendlines)
  cat("\n Elementary components included in the estimation \n"); cat("\n")
  print(x$m)
}

plot.msst <- function(x, time.format="%m-%Y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                        ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL, cex.lab = NULL, cex.axis = NULL, cex.main = NULL,   
                         options = list(type = 'trendlines', 
                                        ncomp.scree = NULL, 
                                        ncomp.erc = rep(1,x$trendlines$D)), ...) {
  if(options$type %!in% c('trendlines', 'screeplot', 'components', 'cpgrams')) {
    print('options must be one of the strings: "trendlines", "components", "cpgrams", or "screeplot".')}
  
 if(options$type == 'trendlines') {
      matplot(x$trendlines$date, x$trendlines$Y,
              main = if(is.null(main)){''} else {main},
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              type = if(is.null(type)){'b'} else {type},
              lwd = if(is.null(lwd)){1} else {lwd},
              xlab = if(is.null(xlab)){'Time'} else {xlab},
              ylab = if(is.null(ylab)){'Singular Spectrum Trendlines'} else {ylab}, 
              ylim = if(is.null(ylim)){c(min(x$trendlines$Y),max(x$trendlines$Y))} else {ylim}, 
              xlim = if(is.null(xlim)){c(min(x$trendlines$dates),max(x$trendlines$dates))} else {xlim}, 
           cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
          cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
          cex.main = if(is.null(cex.main)){1} else {cex.main},
              xaxt = "n")
   if(inherits(x$trendline$dates, "Date")==T){  timelabels<-format(x$trendline$dates,time.format) ;     axis(1,at=x$trendline$dates,labels=timelabels)} else {  axis(1,at=x$trendline$dates)    } 
  }
  
  if(options$type == 'screeplot') {
    if(is.null(options$ncomp.scree)){options$ncomp.scree = (1:3)}
      plot(options$ncomp.scree, log(x$svd$d[options$ncomp.scree]),
           ylab = if(is.null(ylab)){'Eigenvalues (log-scale)'} else {ylab}, 
           xlab = if(is.null(xlab)){'Index'} else {xlab}, 
           lty  = if(is.null(lty)){1} else {lty},
           type = if(is.null(type)){'b'} else {type},
           main = if(is.null(main)){''} else {main},
           pch  = if(is.null(pch)){1} else {pch},
           col = if(is.null(col)){'black'} else {col}, #,
           xlim = if(is.null(xlim)){c(0,max(options$ncomp.scree))} else {xlim},
           ylim = if(is.null(ylim)){c(min(log(x$svd$d[options$ncomp.scree])),
                                      max(log(x$svd$d[options$ncomp.scree])))} else {ylim},
        cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
       cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
       cex.main = if(is.null(cex.main)){1} else {cex.main} )

  }

 if(options$type == 'components') {
      if(is.null(options$ncomp)){options$ncomp = 1 ; print('ncomp is missing, only the first erc is ploted')}
      if(max(options$ncomp) > max(x$m) ){stop('Incompatible number of elementary reconstructed components')}  
      coldim <- options$ncomp; rowdim <- x$trendlines$D     
      par(mfrow = c(rowdim, length(coldim)))
      for(i in 1:rowdim) {
               for(j in coldim) {
            plot(x$trendline$dates, x$erc[[i]][, j],
             main = if(is.null(main)){paste("Serie: ",i," - ERC:", j)} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){''} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c( min(x$erc[[i]][, j]),max(x$erc[[i]][, j]) ) } else {ylim},
             xlim = if(is.null(xlim)){c(min(x$trendlines$dates),max(x$trendlines$dates))} else {xlim}, 
              cex.lab = if(is.null(cex.lab)){0.8} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){0.8} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main} )
                                }  
                           }
 }

if(options$type == "cpgrams") {
  D <- x$trendlines$D
  for(i in 1:D){  cpgram(x$residuals$Y[,i],
                         main = if(is.null(main)){""} else {paste(options$series.names[i])})  }
                              }   
}
