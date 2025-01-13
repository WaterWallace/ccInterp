inversecumtrapz <- function( x, accum, startpos)
{
  # accum is the trapezoidal method of integration
  decum <- c(0, diff(accum))
  recreated_y <- rep(0, length(decum))
  recreated_y[1] <- startpos
  for(i in 1:(length(decum)-1) )
  {
    recreated_y[i+1] <- ( ( decum[i+1] * 2 ) / diff(x)[i] ) - recreated_y[i]
  }
  return(recreated_y)
}
