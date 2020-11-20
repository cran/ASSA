misst <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE)
  UseMethod("misst")

misst.default <- function(y, l = 'automatic', m = 'automatic',
                          vertical = TRUE) {
n <- y$n
D <- y$D
## H. Hassani & R. Mahmoudvand (2013; pp. 68, eq. 25).
if((l == 'automatic') & (vertical == T | vertical == TRUE)) {l <- ceiling((n + 1) / (D + 1))} 
if((l == 'automatic') & (vertical == F | vertical == FALSE)){l <- ceiling(D * (n + 1) / (D + 1))} 
k <- n - l + 1
  
## Run a basic input validation
if (class(y) != "mitsframe")
  stop('Stop: the input y must be a multivariate interval valued object (see ?mitsframe).')
if (is.numeric(l) & (l%%1 != 0 | l <= 0 | l > y$n))
  stop('Stop: l must be a positive integer smaller than number of observations per series.')
if(is.numeric(m) & length(m) != D)
  stop('Stop: The length of vector m must be equal to the number of time series to analyze.')
  
  
## Step 1: 
Y_list = list();
for(d in 1:D){
    Y <- array(NA, c(l, k, 2) ) # l rows, k columns, 2 dimensions
    Y[ , , 1] <- trajectory(y$A[, d], l, k)
    Y[ , , 2] <- trajectory(y$B[, d], l, k) 
    Y_list[[d]] <- Y
  }
  
if(vertical == T | vertical == TRUE) {
  Y <-  array(NA, c(D*l, k,2) )
    for( i in 1:D){   Y[((i-1)*l+1):(i*l) , ,] <- Y_list[[i]] }  } else {
      Y <-  array(NA, c(l, D*k,2) )
      for( i in 1:D){   Y[ ,((i-1)*k+1):(i*k), ] <- Y_list[[i]] }         }
  
## Step 2: SVD
if(vertical == T | vertical == TRUE) {
  S <-  matrix(NA, nrow = D*l, ncol = D*l)
    for( i in 1:D){
      for(j in 1:D){
        S[((i-1)*l+1):(i*l),((j-1)*l+1):(j*l)] <- Var.Est(Y_list[[i]],Y_list[[j]])$outer.product;  
      }
    } } else {
      S <-  matrix(0, nrow = l, ncol = l)
      for( i in 1:D){
        S = S + Var.Est(Y_list[[i]],Y_list[[i]])$outer.product 
      }
} # isSymmetric(S) 
  
SVD <- svd(S)
rk <- qr(SVD$u)$rank # rank.

## Step 3 and 4:
A.tilde = B.tilde = A.Residuals = B.Residuals = matrix(0, ncol = D, nrow = n)
if(is.numeric(m) == T) {
  if(sum(m > rk) > 0) {
    m = replace(m, m > rk, rk)
    warning( "Number of components (see 'm') required to construct trendlines in at least 
                   ^ one dimension was reduced to be compatible with the rank of the trajectory matrix." ) } 
if(vertical == TRUE){
    for(d in 1:D){
      Y_i = YYI_i(Y, SVD, i=m[d], l=D*l, k=k) ;
      tempA           <- dbar(Y_i$Y_iA[((d - 1) * l + 1): (d * l), ], l = l, k = k)
      tempB           <- dbar(Y_i$Y_iB[((d - 1) * l + 1): (d * l), ], l = l, k = k)
      A.tilde[,d]     <- pmin(tempA, tempB)
      B.tilde[,d]     <- pmax(tempA, tempB)
      A.Residuals[,d] <- pmin(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d]) 
      B.Residuals[,d] <- pmax(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d])
    } }   else {
      for(d in 1:D){
        Y_i = YYI_i(Y, SVD, i=m[d], l=l, k=D*k) ; 
          tempA     <- dbar(Y_i$Y_iA[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k) 
          tempB     <- dbar(Y_i$Y_iB[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k)
          A.tilde[,d]     <- pmin(tempA, tempB)
          B.tilde[,d]     <- pmax(tempA, tempB)
          A.Residuals[,d] <- pmin(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d]) 
          B.Residuals[,d] <- pmax(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d])
  } } } else{ 
    m = rep(0,D); #* Storing here the number of ERC on each dimension.
    for(d in 1:D){
      stop.flag = 1
      while(stop.flag != 0 & m[d] < rk) { 
        m[d] = m[d] + 1
        if(vertical == TRUE){
        Y_i = YYI_i(Y, SVD, i=m[d], l=D*l, k=k) ; 
        tempA     <- dbar(Y_i$Y_iA[((d - 1) * l + 1): (d * l), ], l = l, k = k)
        tempB     <- dbar(Y_i$Y_iB[((d - 1) * l + 1): (d * l), ], l = l, k = k)
        A.tilde[,d]     <- pmin(tempA, tempB)
        B.tilde[,d]     <- pmax(tempA, tempB)
        A.Residuals[,d] <- pmin(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d]) 
        B.Residuals[,d] <- pmax(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d])
      }    else {
          Y_i = YYI_i(Y, SVD, i=m[d], l=l, k=D*k) ; 
          tempA     <- dbar(Y_i$Y_iA[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k) 
          tempB     <- dbar(Y_i$Y_iB[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k)
          A.tilde[,d]     <- pmin(tempA, tempB)
          B.tilde[,d]     <- pmax(tempA, tempB)
          A.Residuals[,d] <- pmin(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d]) 
          B.Residuals[,d] <- pmax(y$A[,d]-A.tilde[,d],y$B[,d]-B.tilde[,d])
        } 
  stop.flag <- ecip(as.matrix(cbind(A.Residuals[,d],B.Residuals[,d])));  
} } } # end 'else'

#### ERC:
erc = list()
  for(j in 1:D){
  erc.d <- array(NA, c(n, m[j], 2) ) # n rows, m columns, 2 dimensions
  for(i in 1:m[j]){      
    if(vertical == TRUE){
        Y_i = YYI_i(Y, SVD, i=i, l=D*l, k=k) ;
        ll           <- dbar(Y_i$Y_iA[((d - 1) * l + 1): (d * l), ], l = l, k = k)
        uu           <- dbar(Y_i$Y_iB[((d - 1) * l + 1): (d * l), ], l = l, k = k)
        erc.d[,i,1] <- pmin(ll,uu); erc.d[,i,2] <- pmax(ll,uu)}    else {
        Y_i = YYI_i(Y, SVD, i=i, l=l, k=D*k) ; 
          ll     <- dbar(Y_i$Y_iA[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k) 
          uu     <- dbar(Y_i$Y_iB[ , ((d - 1) * k + 1) : (d * k) ], l = l, k = k)
          erc.d[,i,1] <- pmin(ll,uu); erc.d[,i,2] <- pmax(ll,uu)        } }
  erc[[j]] <- erc.d 
}

# Grouping the results into mtsframe objects to deliver in the output
Y.tilde   <- mitsframe(dates = y$dates, A = A.tilde,     B = B.tilde)
residuals <- mitsframe(dates = y$dates, A = A.Residuals, B = B.Residuals)
  
outputs <- list(trendlines = Y.tilde,
                l = l, 
                m = m, 
                vertical = vertical,
                residuals = residuals,
                svd = SVD,
                erc = erc,
                observations = y,
                call = match.call()) ;  
class(outputs) <- "misst"
return(outputs)
}

print.misst <- function(x, ...) {
  cat("\n Multivariate Singular Spectrum Trendlines for Interval Data:\n ========================================= \n")
  print(x$call)
  cat("\n Interval Trendlines:\n"); cat("\n")
  print(x$trendlines)
  cat("\n Elementary components included in the estimation \n"); cat("\n")
  print(x$m)
}

plot.misst <- function(x, time.format = "%m-%y", col = NULL, lty = NULL,  main = NULL, type = NULL, pch = NULL, lwd = NULL,
                       ylab = NULL,xlab = NULL, ylim = NULL, xlim = NULL, cex.lab = NULL, cex.axis = NULL, cex.main = NULL,
                       tick = FALSE,  options = list(type = 'trendlines' , ncomp = 1:5), ...) {
  if(options$type %!in% c('trendlines', 'screeplots', 'components', 'cpgrams')) {
    print('options must be one of the strings: "trendlines", "components", "cpgrams", or "screeplots".')}
  
  if(options$type == 'trendlines') {
    matplot(x$trendlines$dates, x$trendlines$A,
            main = if(is.null(main)){''} else {main},
            col  = if(is.null(col)){'black'} else {col},
            lty  = if(is.null(lty)){1} else {lty},
            pch  = if(is.null(pch)){1} else {pch},
            lwd  = if(is.null(lwd)){1} else {lwd},
            type = if(is.null(type)){'l'} else {type},
            xlab = if(is.null(xlab)){'Time'} else {xlab},
            ylab = if(is.null(ylab)){'Interval Singular Spectrum Trendline'} else {ylab}, 
            ylim = if(is.null(ylim)){c(min(x$trendlines$A,x$trendlines$B),max(x$trendlines$A,x$trendlines$B))} else {ylim}, 
            xlim = if(is.null(xlim)){c(min(x$trendlines$dates),max(x$trendlines$dates))} else {xlim}, 
            cex.lab = if(is.null(cex.lab)){1} else {cex.lab},
            cex.axis = if(is.null(cex.axis)){1} else {cex.axis},
            cex.main = if(is.null(cex.main)){1} else {cex.main},
            xaxt = "n")
    if(inherits(x$trendlines$dates, "Date")==T){  timelabels<-format(x$trendlines$dates,time.format) ;      
    axis(1,at=x$trendlines$dates,labels=timelabels, tick = tick, cex.axis = if(is.null(cex.axis)){1} else {cex.axis})} else {  axis(1,at=x$trendlines$dates, tick = tick, cex.axis = if(is.null(cex.axis)){1} else {cex.axis})    }
    
    for(i in 1:x$trendlines$D){
      lines(x$trendlines$dates, x$trendlines$A[,i],
            col  = if(is.null(col)){'black'} else {col},
            lty  = if(is.null(lty)){1} else {lty},
            pch  = if(is.null(pch)){1} else {pch},
            lwd  = if(is.null(lwd)){1} else {lwd},
            type = if(is.null(type)){'l'} else {type})
      lines(x$trendlines$dates, x$trendlines$B[,i], 
            col  = if(is.null(col)){'black'} else {col},
            lty  = if(is.null(lty)){1} else {lty},
            pch  = if(is.null(pch)){1} else {pch},
            lwd  = if(is.null(lwd)){1} else {lwd},
            type = if(is.null(type)){'l'} else {type})
      polygon(c(rev(x$trendlines$dates), x$trendlines$dates),
              c(rev(x$trendlines$A[,i]), x$trendlines$B[,i]), 
              col  = if(is.null(col)){'lightgray'} else {col}, border = NA)}
  }
  
  if(options$type == 'screeplots') {
    if(is.null(options$ncomp)){options$ncomp = c(1:3)}
    D <- x$trendlines$D
    for(i in 1:D) {
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
  }
  
  
  if(options$type == 'components') {
    if(is.null(options$ncomp)){options$ncomp = 1 ; print('ncomp is missing, only the first erc is ploted')}
    coldim <- options$ncomp; rowdim <- x$trendlines$D     
    par(mfrow = c(rowdim, length(coldim)))
    for(i in 1:rowdim) {
      for(j in coldim) {
        plot(x$trendline$dates, x$erc[[i]][, j, 1],
             main = if(is.null(main)){paste("Serie: ",i," - ERC:", j)} else {main},
             col  = if(is.null(col)){'black'} else {col},
             lty  = if(is.null(lty)){1} else {lty},
             pch  = if(is.null(pch)){1} else {pch},
             lwd  = if(is.null(lwd)){1} else {lwd},
             type = if(is.null(type)){'l'} else {type},
             xlab = if(is.null(xlab)){''} else {xlab},
             ylab = if(is.null(ylab)){''} else {ylab},
             ylim = if(is.null(ylim)){ c(min(c(x$erc[[i]][, j, 1],x$erc[[1]][, j, 2])),max(c(x$erc[[i]][, j, 1],x$erc[[i]][, j, 2]))) } else {ylim},
             cex.lab = if(is.null(cex.lab)){0.8} else {cex.lab},
             cex.axis = if(is.null(cex.axis)){0.8} else {cex.axis},
             cex.main = if(is.null(cex.main)){1} else {cex.main} )
        lines(x$trendline$dates, x$erc[[i]][, j, 2],  
              col  = if(is.null(col)){'black'} else {col},
              lty  = if(is.null(lty)){1} else {lty},
              pch  = if(is.null(pch)){1} else {pch},
              lwd  = if(is.null(lwd)){1} else {lwd},
              type = if(is.null(type)){'l'} else {type})
        polygon(c(rev(x$trendline$dates), x$trendline$dates),
                c(rev(x$erc[[i]][, j, 1]), x$erc[[i]][, j, 2]), 
                border = NA, col = if(is.null(col)){'lightgray'} else {col}) 
      }  
    }
  }
  
  if(options$type == "cpgrams") {
    D <- x$trendlines$D
    for(i in 1:D){
      ecip(residuos = cbind(x$residuals$A[,i],x$residuals$B[,i]), plot.flag = TRUE)
    }    
  }
}