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

combplot <- function(fit) {
    selection <- rep(0, fit$l)
    selection[fit$sfisher] <- 1
    s <- length(fit$sfisher)
        
    par(mfrow = c(1, 2))
    plot(1, selection[1], type = "h", xlab = "Principal Component",
         ylab = expression(paste("Fisher ", italic(g), " Indicator")),
         lwd = 3, col = "blue", xlim = c(0, fit$l), ylim = c(0, 1))
    for(i in 2:fit$l)
        lines(i, selection[i], type = "h", 
              lwd = 3, col = colorRampPalette(c("blue", "violet",
                                                "red", "orange"))(fit$l)[i])
    plot(fit$erc[, 1], type = "l", lwd = 1, col =
              colorRampPalette(c("blue", "violet", "red",
              "orange"))(fit$l)[fit$sfisher[1]], ylab = "")
    for(i in 2:s)
        lines(fit$erc[, i], type = "l", lwd = 1, col =
              colorRampPalette(c("blue", "violet", "red",
                                 "orange"))(fit$l)[fit$sfisher[i]])
}
