sst <- function(y, l = 'automatic', m = 'automatic')
  UseMethod("sst")

sst.default <- function(y, l = 'automatic', m = 'automatic') {  
    n <- dim(y$Y)[1]
    D <- dim(y$Y)[2]
    Y_ <- matrix(NA, nrow = n, ncol = D) # Estimated trendlines.
    if(is.numeric(m) & length(m) != D) {
        print('Please select the number of components to
be considered for all time series under analysis.')
    }
    
    selection.criteria <- c()
    rank <- c()
    selected.components <- c()
    signals.temp = list()
    eigen.val <- list()
    l.end <- c()
    if(is.numeric(l) & is.numeric(m)) {
        for(i in 1:D) {
            yi <- mtsframe(dates = y$dates, Y = y$Y[, i])
            ti <-  msst(yi, l = l[i], m = m[i])
            Y_[, i] <- ti$trendlines$Y
            selection.criteria[i] <- ti$selection.criteria
            rank[i] <- ti$rank
            l.end[i] <- ti$l
            selected.components[i] <- ti$selected.components
            signals.temp[[i]] <-  ti$erc
            eigen.val[[i]] <-ti$eigen.val
        }
    } else if(is.numeric(l) & !is.numeric(m)) {
        for(i in 1:D) {
            yi <- mtsframe(dates =y$dates, Y = y$Y[, i] )
            ti <-  msst(yi, l = l[i], m = 'automatic')
            Y_[, i] <- ti$trendlines$Y
            selection.criteria[i] <- ti$selection.criteria
            rank[i] <- ti$rank
            l.end[i] <- ti$l
            selected.components[i] <- ti$selected.components
            signals.temp[[i]] <-  ti$erc
            eigen.val[[i]] <-ti$eigen.val
        }
    } else if(!is.numeric(l) & is.numeric(m)) {
        for(i in 1:D) {
            yi <- mtsframe(dates = y$dates, Y = y$Y[, i])
            ti <-  msst(yi, l = "automatic", m = m[i])
            Y_[, i] <- ti$trendlines$Y
            selection.criteria[i] <- ti$selection.criteria
            rank[i] <- ti$rank
            l.end[i] <- ti$l
            selected.components[i] <- ti$selected.components
            signals.temp[[i]] <- ti$erc
            eigen.val[[i]] <-ti$eigen.val
        } 
    } else {
        for(i in 1:D){
            yi <- mtsframe(dates =y$dates, Y = y$Y[, i] )
            ti <-  msst(yi, l = "automatic", m = "automatic")
            Y_[, i] <- ti$trendlines$Y
            selection.criteria[i] <- ti$selection.criteria
            rank[i] <- ti$rank
            l.end[i] <- ti$l
            selected.components[i] <- ti$selected.components
            signals.temp[[i]] <-  ti$erc
            eigen.val[[i]] <-ti$eigen.val
        }    
    }
    
    E_ <- Y_ - y$Y # Matrix with the residuals of each ts
    
    y_ <- mtsframe(dates = y$date, Y = Y_)
    e_ <- mtsframe(dates = y$date, Y = E_)
    signals <- list()
    
    for(i in 1:D)
        signals[[i]] <- do.call(cbind, signals.temp[[i]])
    
    outputs <- list(trendlines = y_, errors = e_, erc = signals,
                    eigen.val = eigen.val, l = l.end,
                    selected.components = selected.components,
                    selection.criteria = selection.criteria,    
                    rank = rank, call = match.call()) ; 
    class(outputs) <- "sst"
    return(outputs)
} 

print.sst <- function(x, ...) {
    cat("\n Singular Spectrum Trendline:\n ============================ \n")
    print(x$call)
    cat("\n Trendlines:\n"); cat("\n")
    print(x$trendlines)
    cat("\n Elementary components included in the estimation \n"); cat("\n")
    print(x$selected.components)
    ## cat("\n Flags indicating the rejection of the null hypothesis (white noise errors): \n"); cat("\n")
    ## print(x$selection.criteria)
}

plot.sst <- function(x, xlab = "Time", ylab = "Singular Spectrum Trendline",
                     options = list(type = "trendlines",
                                    serie = 1, ncomp = 1, series.names = NULL), ...) {
    if(options$type %!in% c("trendlines", "screeplots",
                            "components", "cpgrams")) 
        print("options must be one of the strings:
'trendlines', 'components', 'components', 'cpgrams', or 'screeplots'.") 
    if(options$type=="trendlines") {    
        if(inherits(x$trendlines$dates, "Date")) {      
            matplot(x$trendlines$date, x$trendlines$Y, 
                    ylab = "Singular Spectrum Trendline",
                    type = "l", ylim = c(min(x$trendlines$Y),
                                         max(x$trendlines$Y)), 
                    xaxt = "n", ...)
            timelabels <- format(x$trendlines$date,"%m-%Y")
            axis(1, at = x$trendlines$date, labels = timelabels,
                 tck = 0)                
        } else {
            matplot(x$trendlines$date, x$trendlines$Y, 
                    ylab = "Singular Spectrum Trendline",
                    type = "l", ylim = c(min(x$trendlines$Y),
                                         max(x$trendlines$Y)),...)}
    }
    
  if(options$type == 'screeplots') {
    D <- dim(x$trendlines$Y)[2]
    if(is.list(x$eigen.val)) {
      for(i in 1:D) {
        oldpty <- par(pty = "s")
        on.exit(par(oldpty))
        
        if( is.null(options$series.names) ){
          plot(log(x$eigen.val[[i]][options$ncomp]),
               ylab = 'Eigenvalues (log-scale)', 
               xlab = '', ...) } else {
                 plot(log(x$eigen.val[[i]][options$ncomp]),
                      ylab = 'Eigenvalues (log-scale)', 
                      xlab = '',main=paste(options$series.names[i]), ...)
               }              
        
      } # end for.
    } else {
      oldpty <- par(pty = "s")
      on.exit(par(oldpty))
      if( is.null(options$series.names) ){
        plot(log(x$eigen.val[options$ncomp]),
             ylab = 'Eigenvalues (log-scale)', 
             xlab = '', ...) } else {
               plot(log(x$eigen.val[options$ncomp]),
                    ylab = 'Eigenvalues (log-scale)', 
                    xlab = '', main=paste(options$series.names[i]), ...)
             }
    }
}
    
    if(options$type == "components") {
        if(length(options$ncomp) > dim(x$erc[[options$serie]])[2] |
           max(options$ncomp) > dim(x$erc[[options$serie]])[2]) {
            print(paste("Please choose the number of elementary
components for ploting trendline", options$serie,
" to in the range from 1 to", dim(x$erc[[options$serie]])[2]))  
        } else {            
            coldim <- options$ncomp
            rowdim <- 1            
            par(mar = c(2, 3, 1, 0.5), mfrow = c(rowdim, length(coldim)))
            for(j in coldim) {
                plot(x$erc[[options$serie]][,j], type = "l",
                     ylab = "",
                     main = paste("ERC -",j," Serie", options$serie),
                     cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8)}
        }
    }
 
   if(options$type == "cpgrams") {
    D = dim(x$errors$Y)[2]
    if( is.null(options$series.names) ){
      for(i in 1:D) cpgram(x$errors$Y[,i],main = "") } else {
        for(i in 1:D) cpgram(x$errors$Y[,i],main = paste(options$series.names[i]) )  
      }    
  }
}

