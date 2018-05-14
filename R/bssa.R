##  ========================================================================  ##
##  Miguel de Carvalho                                                        ##
##  Copyright (C) 2018                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

bssa <- function(y, l = 32)
    UseMethod("bssa")

bssa.default <- function(y, l = 32) {
    n <- length(y)
    ## Run a basic input validation
    if (is.ts(y) == FALSE | is.null(dim(y)) == FALSE)
        stop('y must be a univariate time series (ts) object')
    if (l%%1 != 0 | l <= 0 | l > n)
        stop('l must be a positive integer smaller than number of observations')

    k <- n - l + 1
    ## Embedding
    Y <- trajectory(y, l, k)
    ## Singular value decomposition
    SVD <- svd(Y)
    U <- SVD$u
    PC <- matrix(0, nrow = n, ncol = l)
    sfisher <- numeric()
    ## Targeted grouping
    for(i in 1:l) {
        PC[, i] <- dbar(U[, i]%*%t(U[, i])%*%Y, l, k)
        sp <- spec.pgram(ts(PC[, i], frequency = frequency(y),
                            start = start(y)), plot = FALSE, taper = 0.5)
        w <- sp$freq
        spec <- sp$spec
        wstar <- w[which(spec == max(spec))]
        J <- length(w)
        g <- max(spec) / sum(spec)
        aux <- rep(0, J)
        for(j in 1:J)
            aux[j] <- (-1)^(j - 1) * choose(J, j) *
                               max((1 - j * g), 0)^(J - 1)
        pval <- sum(aux)
        if (pval < 0.05 & wstar > 4 / 32 & wstar < 4 / 6) 
            sfisher <- c(sfisher, i)        
    }

    erc <- PC[, sfisher]

    cycle <- as.numeric(rowSums(erc))
    cycle <- ts(cycle, frequency = frequency(y), start = start(y))
    erc <- ts(erc, frequency = frequency(y), start = start(y))
    ## Organize and return outputs
    outputs <- list(cycle = cycle, sfisher = sfisher, erc = erc, l = l, 
                    call = match.call())
    class(outputs) <- "bssa"
    return(outputs)
}

print.bssa <- function(x, ...) {
    cat("\n Singular Spectrum Business Cycle Analysis:\n ==========================================\n ")
    print(x$call)
    cat("\n Business cycle indicator:\n")
    cat("\n")
    print(x$cycle)
    cat("\n Principal components selected by the Fisher g statistic \n")
    cat("\n")
    print(x$sfisher)
}

plot.bssa <- function(x, ylab = "Singular Spectrum Indicator", 
                      lwd = 3, ...) {
    plot(x$cycle, ylab = "Singular Spectrum Indicator", ...)
}

