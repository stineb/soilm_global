stress_exp <- function( x, y0, curve ){
	outstress <- (1-y0) * -exp( -curve * (x - 0.05) ) + 1
	return( outstress )
}

stress_exp_alpha <- function( x, alpha, apar_y0, bpar_y0, apar_curve, bpar_curve ){

	y0 <- apar_y0 + bpar_y0 * alpha
	curve <- apar_curve + bpar_curve * y0

	outstress <- (1 - y0) * -exp( -curve * (x - 0.05) ) + 1

  return( outstress )
}
