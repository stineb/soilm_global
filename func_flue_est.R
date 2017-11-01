## get flue0 : maximum fLUE reduction at low soil moisture
calc_flue0 <- function( alpha ){ 
  ## equation from correct_bysoilm.R
  y <- 0.1 + 0.92 * alpha
  y <- min( y, 1.0 )
  y <- max( y, 0.0 )
  return( y )
}

calc_beta <- function( flue0 ){
  x0 <- 0.125
  x1 <- 0.75
  y <- (flue0 - 1.0) / (x0 - x1)^2
  return( y )
}

calc_flue_est <- function( soilm, beta ){

  x1 <- 0.75

  # if (soilm>x1) {
  #   y <- 1
  # } else {
  #   y <- beta * ( soilm - x1 )^2 + 1.0
  # }
  # y <- min( 1.0, beta * ( soilm - x1 )^2 + 1.0 )
  y <- beta * ( soilm - x1 )^2 + 1.0

  return( y )

}