#### MSSA core functions
trajectory <- function(x, nrow = length(x)%/%2, ncol = length(x)%/%2) {
    # Trajectory matrix 
    Z  <-  matrix(x[1:ncol], nrow = nrow, ncol = ncol, byrow = T)
    for (i in 1:(nrow - 1)) 
	Z[i + 1, ] <- x[-(1:i)][1:ncol]
    Z
}

dbar <- function(Y, l, k) {
    ## Defines R wrapper function to call the FORTRAN function DBAR    
    n <- k + l - 1
    run <- .Fortran('dbar_', Y = as.matrix(Y, nrow = l, ncol = k),
                    l = as.integer(l), k = as.integer(k),
                    answer = numeric(n), PACKAGE = "ASSA")
    return(run$answer)  
}
################################################################
### Signal decomposition when m is a numeric vector       ######
################################################################
## Y = vertical trajectory matrix;
## SVD = svd decomposition of Y%*%t(Y)
#  n = series length. 
#  D = number of series.
#  l = window length.
#  k = n - l + 1
#  m = a vector with the number of components on each dimension.
#  vertical = flag variable; if vertical or horizontal bind is used.
sig.dec1 <- function(Y, SVD, n, D, l, k, m, vertical) {
  e <- SVD$d; U <- SVD$u; rk <- max(m) # eigenvalues and left eigenvectors
  
  # signal.temp a matrix of dimensions "(n,D)" cooresponging to signal q  
  signals <- list()

    if(vertical == T | vertical == TRUE) {
        for(i in 1:rk) {
            v_i <- t(Y)%*%U[, i]/sqrt(e[i])
            S_i <- sqrt(e[i]) * U[, i]%*%t(v_i)
            signals.temp <- matrix(NA, nrow = n, ncol = D)
            for(j in 1:D) {
                temp <- S_i[c((j-1)*l + 1):c(j*l), ]
                signals.temp[, j] <- dbar(temp, l, k)
            }
            signals[[i]] <- signals.temp

        }
    } else {
        for(i in 1:rk) {
            v_i <- t(Y)%*%U[,i] / sqrt(e[i])
            S_i <- sqrt(e[i]) * U[,i]%*%t(v_i)
            signals.temp <- matrix(NA, nrow = n, ncol = D)
            for(j in 1:D) {
                temp <- S_i[  ,c((j-1)*k + 1) : c(j*k)   ]
                signals.temp[ , j] <- dbar(temp, l, k)
            }
            signals[[i]]<-signals.temp 
        }
    } #if else end.
  return(signals) 
}  

##########################################################################
####### An efficient function to extract signals when m = 'automatic'  ###
##########################################################################
sig.dec2 <- function(Y, SVD, n, D, l, k, m, vertical) {
    q <- max(m)
    # signal.temp a matrix of dimensions "(n,D)" cooresponging to signal q  
    signals.temp <- matrix(NA, nrow = n, ncol = D)
  
    v_i <- t(Y)%*%SVD$u[,q] / sqrt(SVD$d[q])
    S_i <- sqrt(SVD$d[q]) * SVD$u[,q]%*%t(v_i)
  
    if(vertical == T | vertical == TRUE) {
        for(j in 1:D) {
          temp <- S_i[c((j-1)*l + 1):c(j*l), ]
          signals.temp[ , j] <- dbar(temp, l, k)
        }           
    } else {
        for(j in 1:D) {
            temp <- S_i[, c((j-1)*k + 1):c(j*k)]
            signals.temp[, j] <- dbar(temp, l, k)} 
    }      
    return(signals.temp) 
}  
##########################################################################
 '%!in%' <- function(x,y) {!('%in%'(x,y))}
##########################################################################
## Projections onto the simplex:
## Input : the vector p      in R*n
## Output: the vector p.simp in S*n
simp.proj <- function(p) {
    p <- as.vector(p)
    p.sort <- sort(p, decreasing = TRUE);        
    n <- length(p.sort)
    j <- 1
    while(p.sort[j] + (1 - sum(p.sort[1:j])) / j > 0 & j <= n)
        j <- j + 1
    j <- j - 1
    lambda <- (1 - sum(p.sort[1:j])) / j
    p.simp <- rep(0, n)
    for(i in 1:n)
        p.simp[i] <- max(p[i] + lambda, 0)
    return(p.simp)
}

###################################################################################
###  This functions is a remake of cpgram already available on r-package stats.  ##
###################################################################################
cpgram2 <- function(ts, taper = 0.1, plot = FALSE, 
                    main = paste("Series: ", deparse(substitute(ts))), 
                    ci.col = "blue") {  
    x <- as.vector(ts)
    x <- x[!is.na(x)]
    x <- spec.taper(scale(x, TRUE, FALSE), p = taper)
    y <- Mod(fft(x))^2 / length(x)
    y[1L] <- 0
    n <- length(x)
    x <- (0:(n/2)) * frequency(ts)/n
    if (length(x)%%2 == 0) {
        n <- length(x) - 1
        y <- y[1L:n]
        x <- x[1L:n]
    }   else {
        y <- y[seq_along(x)]
    }
    xm <- frequency(ts) / 2
    mp <- length(x) - 1
    ## crit <- 1.628/(sqrt(mp) + 0.12 + 0.11/sqrt(mp)) # p-val = 0.01
    crit <- 1.358/(sqrt(mp) + 0.12 + 0.11 / sqrt(mp)) # p-val = 0.05
  
    ## Our modifications to cpgram starts here:
    ## D_mp is our KS statistic:
    D_mp <- max(abs(cumsum(y) / sum(y) - seq(0, 1, length.out = mp + 1)))
    reject <- as.numeric(D_mp >= crit)
    ## Compute area under the cpgram:
    area.under <- sum(diff(x)[1]*cumsum(y)/sum(y))
    area.crit <- as.numeric(area.under >= xm/2)
    ## out = 1 then reject the null (erros are white noise). 
    if(plot == TRUE) {
        oldpty <- par(pty = "s")
        on.exit(par(oldpty))
        plot(x, cumsum(y)/sum(y), type = "s", xlim = c(0, xm), 
             ylim = c(0,1), xaxs = "i", 
             yaxs = "i", xlab = "frequency", ylab = "")
        lines(c(0, xm * (1 - crit)), c(crit, 1), col = ci.col, lty = 2)
        lines(c(xm * crit, xm), c(0, 1 - crit), col = ci.col, lty = 2)
        title(main = main) 
    }
    if(reject == 1 & area.crit > xm/2) {
        out <- 1
    } else {
        out <- 0
    } 
    return(out = out)
    # and ends here.
}
