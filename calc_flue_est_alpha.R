calc_flue_est_alpha <- function( x, alpha, apar, bpar, cpar, dpar ){
  flue0 <- apar + bpar * alpha
  beta  <- (flue0 - 1.0) / (cpar - dpar)^2
  flue_est <- beta * (x - dpar)^2 + 1
  flue_est <- ifelse( x>dpar, 1, ifelse( x<0, 0, flue_est ) )
  return( flue_est )
}
