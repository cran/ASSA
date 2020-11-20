###############################################################
#### Core functions in ASSA package ###########################
###############################################################

trajectory <- function(x, l = length(x)%/%2, k = length(x)%/%2){
  # Trajectory matrix 
  Z  <-  matrix(x[1:k], nrow = l, ncol = k, byrow = T)
  for (i in 1:(l - 1)) 
    Z[i + 1, ] <- x[-(1:i)][1:k]
  return(Z)
}

Y_i <- function(Y, SVD, i) {
  v_i <- t(Y)%*%SVD$u[,i] / sqrt(SVD$d[i])
  y_i <- sqrt(SVD$d[i])*SVD$u[,i]%*%t(v_i)
  return(y_i) }  # Output: i-th matrix corresponding to the trajectory matrix (Y) decomposition.

YYI_i <- function(Y, SVD, i, l, k) {
  U = matrix(0,l,l)
  for(j in 1:i){U = U + SVD$u[,j]%*%t(SVD$u[,j])}

  ret <- .Fortran('yyi_', YA = as.matrix(Y[ , , 1]), YB = as.matrix(Y[ , , 2]),
                  U = as.matrix(U), l = as.integer(l), k = as.integer(k),    
                  ansA = matrix(0,nrow = l, ncol = k) , 
                  ansB = matrix(0,nrow = l, ncol = k), PACKAGE = "ASSA")
  
    res <- .Fortran('yyi_', YA = as.matrix(Y[ , , 1]), YB = as.matrix(Y[ , , 2]),
                  U = as.matrix(diag(l)-U), l = as.integer(l), k = as.integer(k),    
                  ansA = matrix(0,nrow = l, ncol = k) , 
                  ansB = matrix(0,nrow = l, ncol = k), PACKAGE = "ASSA")
  
  Y_i = list(Y_iA = ret$ansA, Y_iB = ret$ansB, E_A =res$ansA, E_B =res$ansB)
  return(Y_i) }  # Output: list with matrices Yi.

YYI_i.erc <- function(Y, SVD, i, l, k) {
  U = SVD$u[,i]%*%t(SVD$u[,i])
  
  ret <- .Fortran('yyi_', YA = as.matrix(Y[ , , 1]), YB = as.matrix(Y[ , , 2]),
                  U = as.matrix(U), l = as.integer(l), k = as.integer(k),    
                  ansA = matrix(0,nrow = l, ncol = k) , 
                  ansB = matrix(0,nrow = l, ncol = k), PACKAGE = "ASSA")
  
  Y_i = list(Y_iA = ret$ansA, Y_iB = ret$ansB)
  return(Y_i) }  # 

dbar <- function(Y, l, k) {
  ## Defines R wrapper function to call the FORTRAN function DBAR    
  n <- k + l - 1
  run <- .Fortran('dbar_', Y = as.matrix(Y, nrow = l, ncol = k),
                  l = as.integer(l), k = as.integer(k),
                  answer = numeric(n), PACKAGE = "ASSA")
  return(run$answer)  
}

autocov <- function(residuos) {
  ## Defines R wrapper function to call the FORTRAN function AUTOCOV
  n = dim(residuos)[1]
  residuos = scale(residuos, center = TRUE, scale = FALSE)
  run <- .Fortran('autocov_', e = as.matrix(residuos, nrow = n , ncol = 2),
                  n = as.integer(n), answer = numeric(n), PACKAGE = "ASSA")
  return(run$answer)  
}

## Simplex projections:
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


Var.Est <-function(A,B){
  try(if(is.array(A)!=T) stop('Y must be an array'))
  try(if(is.array(B)!=T) stop('X must be an array'))
  try(if(dim(A)[3]!=2) stop('Not suitable dimensions for the array Y, please check.'))
  try(if(dim(B)[3]!=2) stop('Not suitable dimensions for the array X, please check.'))
  l=dim(A)[1];k=dim(A)[2]
  
  outer   <- 2*A[,,1]%*%t(B[,,1]) + A[,,1]%*%t(B[,,2]) + A[,,2]%*%t(B[,,1]) + 2*A[,,2]%*%t(B[,,2])
#  var    <- (2*scale(A[,,1],scale=FALSE)%*%t(scale(B[,,1],scale=FALSE)) + scale(A[,,1],scale=FALSE)%*%t(scale(B[,,2],scale=FALSE)) + 
#              scale(A[,,2],scale=FALSE)%*%t(scale(B[,,1],scale=FALSE)) + 2*scale(A[,,2],scale=FALSE)%*%t(scale(B[,,2],scale=FALSE)))/(6*k)
  output  <- list(outer.product = outer) #
  return(output)
}



## crit <- 1.628/(sqrt(mp) + 0.12 + 0.11/sqrt(mp)) # p-val = 0.01
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
  crit <- 1.358/(sqrt(mp) + 0.12 + 0.11 / sqrt(mp)) # p-val = 0.05
  
  ## Our modifications to cpgram starts here:
  ## D_mp is our KS statistic:
  D_mp <- max(abs(cumsum(y) / sum(y) - seq(0, 1, length.out = mp + 1)))
  reject <- as.numeric(D_mp >= crit) # Rejection flag.
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
  if(reject == 1 & area.crit > xm/2) { out <- 1  } else { out <- 0  } 
  return(out = out)
}


#### ECIP estimation:
ecip <- function(residuos, plot.flag = FALSE){
  mp = dim(residuos)[1]
  omega = unique( (2*pi/mp)*( ceiling( (1:(mp-1))/2 ) ) )
  hat.gamma <- c(autocov(residuos))
  
  M = mp - 1
  I = c() # Periodogram estimation.
  for(i in 1:length(omega)){
    I[i] = hat.gamma[1] + 2*(sum( cos(omega[i]*(1:M))*hat.gamma[2:(M+1)] ))} # BD Eq in Spectral Analysis of Signals (P. Stoica)
  I = c(mp*(mean(colMeans(residuos)))^2 ,I)                                  
  c = sum(abs(I)) # 
  empirical.cdf <- cumsum(abs(I))/c  

  unif.ref <- cumsum(rep(1/length(I),length(I)))
  crit <- 1.36*(ceiling((mp-1)/2) - 1)^{-1/2} # p-val = 0.05
  # crit <- 1.63*(ceiling((mp-1)/2) - 1)^{-1/2} # p-val = 0.01
  flag1 = sum(empirical.cdf > unif.ref + crit)   ### ">0" when the null (withe noise) is rejected.
  flag2 = sum(empirical.cdf < unif.ref - crit)   ### ">0" when the null (withe noise) is rejected.
  
  if(plot.flag == TRUE){
    plot(empirical.cdf, type = 'l', xlab = '', ylab = '',ylim = c(0,1))
    points(unif.ref + crit, type = 'l', col = 'blue', lty = 2)
    points(unif.ref - crit, type = 'l', col = 'blue', lty = 2)
  }
  
  return(max(flag1,flag2))
}

###########################################
### Other internal functions             ##
###########################################
###########################################
'%!in%' <- function(x,y) {!('%in%'(x,y))}
'%notin%' <- Negate('%in%')


########################################### END.
