msstc <- function(y, l = 'automatic', m = 'automatic', vertical = TRUE)
    UseMethod("msstc")

msstc.default <- function(y, l = 'automatic' , m = 'automatic',
                          vertical = TRUE) {
  msst.output <- msst(y = y, l = l, m = m, vertical = vertical)    
  n <- dim(msst.output$trendlines$Y)[1]
  D <- dim(msst.output$trendlines$Y)[2]
  proj.trendline <- matrix(NA, ncol = D, nrow = n)

  for(i in 1:n) 
    proj.trendline[i, ] <- simp.proj(msst.output$trendlines$Y[i,])
                 
  # Grouping the results into mts objects to deliver in the output.
  y1_ <- mtsframe(dates = msst.output$trendlines$date, 
                  Y = proj.trendline)
  
  outputs <- list(trendlines = y1_, errors = msst.output$errors, 
                  erc = msst.output$erc, eigen.val = msst.output$eigen.val,
                  l = l,  selected.components = msst.output$selected.components, 
                  rejection.criteria = msst.output$rejection.criteria,
                  rank = msst.output$rank, call = match.call())  
  class(outputs) <- "msst"
  return(outputs)
}
