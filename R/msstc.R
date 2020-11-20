msstc <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE)
  UseMethod("msstc")

msstc.default <- function(y, l = 'automatic' , m = 'automatic',
                          vertical = TRUE) {
  msst.output <- msst(y = y, l = l, m = m, vertical = vertical)    
  n <- msst.output$trendlines$n
  D <- msst.output$trendlines$D
  proj.trendline <- matrix(NA, nrow = n, ncol = D)
  
  for(i in 1:n) {  proj.trendline[i, ] <- c(simp.proj(msst.output$trendlines$Y[i,])) }
  
  # Grouping the results into mts objects to deliver in the output.
  Y.tilde   <- mtsframe(dates = y$date, Y = proj.trendline)
  Residuals <- mtsframe(dates = y$date, Y = proj.trendline - y$Y)
  
  
  
  outputs <- list(trendlines = Y.tilde, 
                  l = msst.output$l, 
                  m = msst.output$m,
                  vertical = msst.output$vertical, 
                  residuals = Residuals,
                  svd = msst.output$svd,
                  erc = msst.output$erc,
                  observations = msst.output$y,
                  call = match.call()) ;   
  class(outputs) <- "msstc"
  return(outputs)
}

print.msstc <- function(x, ...) {
  cat("\n Multivariate Singular Spectrum Trendlines:\n ========================================= \n")
  print(x$call)
  cat("\n Trendlines:\n"); cat("\n")
  print(x$trendlines)
  cat("\n Elementary components included in the estimation \n"); cat("\n")
  print(x$m)
}

plot.msstc <- function(x, time.format="%m-%Y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
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
