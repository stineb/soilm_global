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

stress_exp_alpha_grasstree <- function( x, alpha, apar_y0, bpar_y0, apar_curve, bpar_curve, classid=NA ){

  if (classid %in% c("GRA", "CSH") ){
    y0 <- apar_y0[2] + bpar_y0[2] * alpha
  } else {
    y0 <- apar_y0[1] + bpar_y0[1] * alpha
  }

  curve <- apar_curve + bpar_curve * y0
  outstress <- (1 - y0) * -exp( -curve * (x - 0.05) ) + 1
  return( outstress )
}