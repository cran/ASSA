predict <- function(fitted.model, p = 1 )
  UseMethod("predict")

predict.default <- function(fitted.model, p = 1) {
  
  if (missing(fitted.model))
    stop('Stop: Please include a model as an imput to produce the forecasting')
  if ( p%%1 != 0 | p <= 0  )
    stop('Stop: The number of periods ahead to forecast must be a positive integer.')
  
  # Parameters:
  
  
  if(class(fitted.model)%in% c('sst')) {
    # the svd's on each dimension are hosted in a list:
    n <- fitted.model$trendline$n
    D <- 1
    m <- fitted.model$m
    l <- fitted.model$l # this is a single number in case of mssa mssac
    y.tilde <- fitted.model$trendline$y
    U   <- as.matrix(fitted.model$svd$u[,1:m])
    a<-rep(NA,l-1)
    
    pi<-U[l,] ; Pi<-as.matrix(U[-l,])
    f <- rep(NA,p); a <-rep(NA,l-1);
    
    v2<-sum(pi^2); aux<-matrix(0,l-1,m)
    
    for(i in 1:m){ aux[,i] <- pi[i] * Pi[,i] }
    a[(l-1):1] <- rowSums(aux) / (1-v2)
    
    for(i in 1:p){
      nn   <- length(y.tilde)
      f[i] <- a %*% y.tilde[nn:(nn-l+2)]
      y.tilde  <- c(y.tilde,f[i]) }
    
    Forecast  <- f; R <- a
    outputs <- list(forecasts = Forecast, coefficients = R)
    class(outputs) <- "predict"
    return(outputs)
  }
  
  
  if(class(fitted.model)%in% c('msst','msstc')) {
    n <- fitted.model$trendlines$n; D <- fitted.model$trendlines$D
    m <- fitted.model$m; l <- fitted.model$l # this is a single number in case of mssa mssac
    r <- max(m)
    Forecast <- matrix(NA, nrow = p, ncol=D)
    ## Creating R matrix according to staking strategy
    if(fitted.model$vertical==TRUE){
      W  <- matrix(0, nrow = D      , ncol = r) # W matrix
      Q <- matrix(0, nrow = (l-1)*D, ncol = r)  # Q matrix (U^\nabbla,M^T) in Hassani).
      Zh <- matrix(0, nrow = (l-1)*D, ncol = 1) # Zh must be row vector (dim: (l-1)D) *Notacion incosistente en Hassani*
      for(j in 1:D){ 
        U   <- as.matrix(fitted.model$svd$u[ ((j-1)*l + 1):(j*l) , (1:r)]) 
        W[j,(1:r)] <- U[l,(1:r)]
        Q[((j-1)*(l-1)+1):(j*(l-1)) ,1:r]  <- U[-l,(1:r)] 
        Zh[((j-1)*(l-1)+1):(j*(l-1)) ,  1 ] <- fitted.model$trendlines$Y[(n-l+2):n, j]
      }
      R <- solve(diag(D) - W%*%t(W))%*%W%*%t(Q) 
      # Recursive forecasting:
      h <- 1; 
      while(h<=p){
        Forecast[h , ] <- t(R%*%Zh)
        # Updates:
        for(j in 1:D){ Zh[((j-1)*(l-1)+1):(j*(l-1)),] <- c(fitted.model$trendlines$Y[,j],Forecast[1:h ,j ])[(n+h-l+2):(n+h)] }
        h   <- h + 1
      } } else {
        # Horizontal staking:
        R  <- matrix(0, nrow = l-1, ncol = 1)
        Zh <- matrix(0, nrow = D, ncol = (l-1))
        U   <- as.matrix(fitted.model$svd$u[ , 1:r])
        v2<-sum(U[l, ]^2);
        for(j in 1:r){R = R + U[l,j]*U[-l,j] }
        R = R/(1-v2)
        for(j in 1:D){ Zh[ j ,  ] <- fitted.model$trendlines$Y[(n-l+2):n, j] }
        # Recursive forecasting:
        h <- 1;
        while(h<=p){
          Forecast[h , ] <- Zh%*%R
          # Updates:
          for(j in 1:D){ Zh[ j ,  ] <- c(fitted.model$trendlines$Y[,j],Forecast[1:h , j ])[(n+h-l+2):(n+h)] }
          h   <- h + 1 } # end 'while'  
      } # end 'else'
    
    
    if(class(fitted.model)%in% c('msstc')) { for(i in 1:p) {  Forecast[i, ] <- c(simp.proj(Forecast[i,])) } }
    
    outputs <- list(forecasts = Forecast, coefficients = R)
    class(outputs) <- "predict"
    return(outputs)
  } # End msst and msstc.  
  
  
  if(class(fitted.model)%in%c('isst') ) {          
    y_ja <- fitted.model$trendline$a
    y_jb <- fitted.model$trendline$b
    m <- fitted.model$m; 
    l <- fitted.model$l;
    U   <-  as.matrix(fitted.model$svd$u[,1:m]) 
    Forecast = matrix(NA, ncol = 2, nrow = p)
    pi<-U[l,] ; Pi<-as.matrix(U[-l,]);
    fa = fb <- rep(NA,p); a<-rep(NA,l-1)
    
    v2<-sum(pi^2); aux<-matrix(0,l-1,m)
    
    for(i in 1:m){ aux[,i] <- pi[i] * Pi[,i] }
    a[(l-1):1] <- rowSums(aux) / (1-v2)
    
    for(i in 1:p){
      nn   <- length(y_ja)
      fa[i] <- a %*% y_ja[nn:(nn-l+2)]
      y_ja  <- c(y_ja,fa[i])
      fb[i] <- a %*% y_jb[nn:(nn-l+2)]
      y_jb  <- c(y_jb,fb[i])
    }
    
    Forecast[ , 1] <- pmin(fa,fb); Forecast[ , 2] <- pmax(fa,fb); R <- a
    outputs <- list(forecasts = Forecast, coefficients = R)
    class(outputs) <- "predict"
    return(outputs)
  } 
  
  if(class(fitted.model)%in% c('misst') ) {    
    n <- fitted.model$trendlines$n
    D <- fitted.model$trendlines$D
    m <- fitted.model$m
    l <- fitted.model$l # this is a single number in case of mssa mssac
    r <- max(m)
    Forecast <- array(NA, c(p, D, 2))
    
    if(fitted.model$vertical==TRUE){
      W  <- matrix(0, nrow = D     , ncol = r) 
      Q <- matrix(0, nrow = (l-1)*D, ncol = r) 
      ZhA <- ZhB  <- matrix(0, nrow = (l-1)*D, ncol = 1) 
      for(j in 1:D){ 
        U   <- fitted.model$svd$u[((j-1)*l+1):(j*l) , (1:r)] 
        W[j,1:r] <- U[l,1:r]
        Q[((j-1)*(l-1)+1):(j*(l-1)) ,1:r] <- U[-l,(1:r)] 
        ZhA[((j-1)*(l-1)+1):(j*(l-1)) ,   1   ] <- fitted.model$trendlines$A[(n-l+2):n, j]
        ZhB[((j-1)*(l-1)+1):(j*(l-1)) ,   1   ] <- fitted.model$trendlines$B[(n-l+2):n, j]
      }
      
      R <- solve(diag(D) - W%*%t(W))%*%W%*%t(Q) 
      # Recursive forecasting:
      h <- 1;
      while(h<=p){
        Forecast[h , ,1] <- t(R%*%ZhA) ;   Forecast[h , ,2] <- t(R%*%ZhB)
        # Updates:
        for(j in 1:D){ ZhA[((j-1)*(l-1)+1):(j*(l-1)),] <- c(fitted.model$trendlines$A[,j],Forecast[1:h , j, 1])[(n+h-l+2):(n+h)];
                       ZhB[((j-1)*(l-1)+1):(j*(l-1)),] <- c(fitted.model$trendlines$B[,j],Forecast[1:h , j, 2])[(n+h-l+2):(n+h)]}
        h   <- h + 1
      } } else {     
      # Horizontal staking:
      R <- matrix(0, nrow = (l-1), ncol = 1)  
      ZhA = ZhB  <- matrix(0, nrow = D, ncol = (l-1))
      U   <- as.matrix(fitted.model$svd$u[ , 1:r])
      v2<-sum(U[l, ]^2);
      for(j in 1:r){R = R + U[l,j]*U[-l,j] }
      R = R/(1-v2)
      for(j in 1:D){ ZhA[ j , ] <- fitted.model$trendlines$A[(n-l+2):n, j]; ZhB[ j , ] <- fitted.model$trendlines$B[(n-l+2):n, j]  }
      # Recursive forecasting:
      h <- 1;
      while(h<=p){
        Forecast[h , , 1] <- ZhA%*%R;  Forecast[h , , 2] <- ZhB%*%R;
        # Updates:
        for(j in 1:D){ ZhA[ j ,  ] <- c(fitted.model$trendlines$A[,j],Forecast[1:h ,j ,1 ])[(n+h-l+2):(n+h)];
                       ZhB[ j ,  ] <- c(fitted.model$trendlines$B[,j],Forecast[1:h ,j ,2 ])[(n+h-l+2):(n+h)]}
        h   <- h + 1 } # end 'while'  
  } # end 'else'
    
  outputs <- list(forecasts = Forecast, coefficients = R)
  class(outputs) <- "predict"
  return(outputs)
  }
}
