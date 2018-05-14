msst <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE)
    UseMethod("msst")

msst.default <- function(y, l = 'automatic', m = 'automatic',
                         vertical = TRUE) {
    n <- y$n
    D <- y$D
    ## Cf H. Hassani & R. Mahmoudvand (2013; pp. 68, eq. 25).
    if((l == 'automatic') & (vertical == T | vertical == TRUE))
        l <- ceiling((n + 1) / (D + 1)) 
    if((l == 'automatic') & (vertical == F | vertical == FALSE))
        l <- ceiling(D * (n + 1) / (D + 1)) 
    k <- n - l + 1
    
    ## Run a basic input validation
    if (class(y) != "mtsframe")
        stop('Stop: the input y must be an mtsframe object.')
    if (is.numeric(l) & (l%%1 != 0 | l <= 0 | l > y$n))
        stop('Stop: l must be a positive integer smaller 
             than number of observations per series.')
    if(is.numeric(m) & length(m) != D)
        stop('Stop: The length of vector m must be equal
to the number of time series.')

    ## Embedding
    if(vertical == T | vertical == TRUE) {
    ## Vertical embedding
        Y <- matrix(NA, nrow = D * l, ncol = k)
        for(d in 1:D)
            Y[((d - 1) * l + 1): (d * l), ] <- trajectory(y$Y[, d], l, k)
        d <- qr(Y)$rank # D = number of series, d = rank of Autocorr matrix
    } else {
    ## Horizontal embedding 
        Y <- matrix(NA, nrow =  l, ncol = D * k)
        for(d in 1:D)
            Y[, ((d - 1) * k + 1) : (d * k) ] <- trajectory(y$Y[, d], l, k)
        d <- qr(Y)$rank  
    }

    ## SVD
    SVD.S <- svd(Y%*%t(Y))
  
    ## Grouping:  
    ## Input validation and consistensy checks:
    if(is.numeric(m) == T) {
        if(sum(m > d) > 0) {
            replace(m, m > d, d)
            warning("Warnings about 'm': the number of components required to 
              construct the estimator in at least 
              one of the series was reduced to be compatible with the
              rank of the trajectory matrix." )
        }    
        if(sum(m < 1) > 0) {
            replace(m, m < 1, 1) 
            warning("Warnings about'm': The number of some components 
               where forced to be at least 1 
              (please, choose only positive integers for 'm' next time)" )
        }
    } 
    
    if(is.numeric(m)) {
        signals <- sig.dec1(Y = Y, SVD = SVD.S, n = n, D = D, l = l,
                            k = k,m = m, vertical = vertical)
        Y_ <- matrix(NA,nrow = n, ncol = D) # Estimated trendlines
  
        for(i in 1:D) {
            if(m[i] > 1) {
                Y_[, i] <- Reduce("+", signals[1:m[i]])[, i] 
            } else {
                Y_[, i] <-signals[[1]][, i]
            }
        }
        
        E_ <- Y_ - y$Y # Matrix with the errors of each ts
        
        # Computing flags for H0 in the errors (on each dimension)
        stop.flag <- rep(NA,D)
        for(i in 1:D){
            stop.flag[i] <- 1-cpgram2(E_[,i])
        }        
    } else {

        stop.flag <- rep(0, D)
        mm <- rep(1,D)
        Y_ <- matrix(NA, nrow = n, ncol = D) # Estimated trendlines.s
        signals = list()
        signal.extraction.counter <- 1

      while(sum(stop.flag) < D) {  
            flag.ind  <- which(stop.flag == 0) # Used later for the trendline        
            signals[[signal.extraction.counter]]  <- sig.dec2(Y = Y, SVD = SVD.S,
                                                              n = n, D = D, l = l,
                                                              k = k, m = mm,
                                                              vertical = vertical)            
            # Producing signals for selected ts
            for(i in flag.ind) {
                if(mm[i] > 1) {
                    Y_[, i] <- Reduce("+", signals[1:mm[i]])[, i] 
                } else {
                    Y_[, i] <-signals[[1]][, i]
                }
            }
        
            E_ <- Y_ - y$Y # Matrix with the errors of each ts
        
           # Updating stoping criterions
            for(i in flag.ind) {
                if(cpgram2(E_[, i]) == 1) {
                    mm[i] <- mm[i] + 1 
                } else {
                    stop.flag[i] <- 1
                }
            }
            
            signal.extraction.counter <- signal.extraction.counter + 1        
      } # while-end
    } 
    
    # Grouping the results into mtsframe objects to deliver in the output
    y_ <- mtsframe(dates = y$date, Y = Y_)
    e_ <- mtsframe(dates = y$date, Y = E_)
    if(is.numeric(m)) {
        selected.comp <- m
    } else {
        selected.comp <- mm
    }
    
    outputs <- list(trendlines = y_, errors = e_, erc = signals, 
                    eigen.val = SVD.S$d, l = l, 
                    selected.components = selected.comp,
                    selection.criteria = 1 - stop.flag,    
                    rank = d, call = match.call())  
    class(outputs) <- "msst"
    return(outputs)
}

print.msst <- function(x, ...) {
    cat("\n Multivariate Singular Spectrum Trendlines:\n ========================================= \n")
    print(x$call)
    cat("\n Trendlines:\n"); cat("\n")
    print(x$trendlines)
    cat("\n Elementary components included in the estimation \n"); cat("\n")
    print(x$selected.components)
    cat("\n Flags indicating the rejection of the
null hypothesis (white noise errors): \n"); cat("\n")
    print(x$selection.criteria)
}

plot.msst <- function(x, xlab = "Time", ylab = "Singular Spectrum Trendline",
                      options = list(type = 'trendlines' , ncomp = 1:5, series.names =NULL ), ...) {
    if(options$type %!in% c('trendlines', 'screeplots', 'components', 'cpgrams')) {
        print('options must be one of the strings:
"trendlines", "components", "components", "cpgrams", or "screeplots".'
)}
    
    if(options$type == 'trendlines') {
        if(inherits(x$trendlines$dates, "Date")) {    
            matplot(x$trendlines$date, x$trendlines$Y, xlab = xlab,  
                    ylab = "Singular Spectrum Trendline",
                    type = "l", ylim = c(min(x$trendlines$Y),
                                         max(x$trendlines$Y)), 
                    xaxt = "n", ...)
            timelabels <- format(x$trendlines$date, "%m-%Y")
            axis(1, at = x$trendlines$date, labels = timelabels, tck = 0) 
        } else {
            matplot(x$trendlines$date, x$trendlines$Y, xlab = xlab, 
                    ylab = "Singular Spectrum Trendline",
                    type = "l", ylim = c(min(x$trendlines$Y),
                                         max(x$trendlines$Y)), ...)
        }
    }
    
    if(options$type == 'screeplots') {
        D <- dim(x$trendlines$Y)[2]
        if(is.list(x$eigen.val)) {
            for(i in 1:D) {
                oldpty <- par(pty = "s")
                on.exit(par(oldpty))
                
                if(is.null(options$series.names) ){
                    plot(log(x$eigen.val[[i]][options$ncomp]),
                         ylab = 'Eigenvalues (log-scale)', 
                         xlab = '', ...)
                } else {
                    plot(log(x$eigen.val[[i]][options$ncomp]),
                         ylab = 'Eigenvalues (log-scale)', 
                         xlab = '', main = paste(options$series.names[i]), ...)
                }                              
            } # end for.
        } else {
            oldpty <- par(pty = "s")
            on.exit(par(oldpty))
            if(is.null(options$series.names)) {
                plot(log(x$eigen.val[options$ncomp]),
                     ylab = 'Eigenvalues (log-scale)', 
                     xlab = '', ...)
            } else {
                plot(log(x$eigen.val[options$ncomp]),
                     ylab = 'Eigenvalues (log-scale)', 
                     xlab = '', main=paste(options$series.names[i]), ...)
            }
        }
    }
    
    if(options$type == 'components') {
        if(length(options$ncomp) > length(x$erc) |
           max(options$ncomp) > length(x$erc)) {
            print( paste('Please choose the elementary
             components to be plot in the range from 1 to', dim(x$erc)[3])) 
       } else {
           coldim <- options$ncomp
           rowdim <- dim(x$trendlines$Y)[2]     
           par(mar = c(2, 3, 1, 0.5), mfrow = c(rowdim, length(coldim)))
           for(i in 1:rowdim) {
               for(j in coldim) {
                   plot(x$erc[[j]][, i],type = "l", ylab = "",
                        main=paste("ERC -", j," Serie", i), cex.lab = 0.8,
                        cex.axis = 0.8, cex.main = 0.8)
               }
           }
       }
    }

    if(options$type == "cpgrams") {
        D <- dim(x$errors$Y)[2]
        if(is.null(options$series.names)){
            for(i in 1:D)
                cpgram(x$errors$Y[, i], main = "")
        } else {
            for(i in 1:D)
                cpgram(x$errors$Y[,i],
                       main = paste(options$series.names[i]))
        }    
    }
}
