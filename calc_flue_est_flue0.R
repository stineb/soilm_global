calc_flue_est_flue0 <- function( x, flue0, cpar, dpar ){
  beta  <- (flue0 - 1.0) / (cpar - dpar)^2
  flue_est <- beta * (x - dpar)^2 + 1
  return( flue_est )
}