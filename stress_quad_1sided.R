stress_quad_1sided <- function( x, x0, beta ){
  outstress <- 1.0 - beta * ( x - x0 ) ^ 2
  outstress <- ifelse( x>x0, 1, outstress )
  return( outstress )
}

stress_quad_1sided_alpha <- function( x, alpha, x0, apar, bpar ){
	y0 <- apar + bpar * alpha
	beta <- (1 - y0) / x0^2
  outstress <- 1.0 - beta * ( x - x0 ) ^ 2
  outstress <- ifelse( x>x0, 1, outstress )
  return( outstress )
}
