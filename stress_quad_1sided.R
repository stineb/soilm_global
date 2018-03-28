stress_quad_1sided <- function( x, x0, beta ){
  outstress <- 1.0 - beta * ( x - x0 ) ^ 2
  outstress <- ifelse( x>x0, 1, ifelse( x<0, 0, outstress ) )
  return( outstress )
}

stress_quad_1sided_alpha <- function( x, alpha, x0, apar, bpar ){

  y0 <- apar + bpar * alpha
  beta <- (1 - y0) / x0^2
  outstress <- 1.0 - beta * ( x - x0 ) ^ 2    
  outstress <- ifelse( x>x0, 1, ifelse( x<0, 0, outstress ) )
  outstress <- ifelse( outstress>1, 1, ifelse( outstress<0, 0, outstress ) )

  return( outstress )

}


stress_quad_1sided_alpha_grasstree <- function( x, alpha, x0, apar, bpar, classid=NA ){

  ## the first element is for non-grasses
  if (classid %in% c("GRA", "CSH") ){
    y0 <- apar[2] + bpar[2] * alpha
  } else {
    y0 <- apar[1] + bpar[1] * alpha
  }

  beta <- (1 - y0) / x0^2
  outstress <- 1.0 - beta * ( x - x0 ) ^ 2    
  outstress <- ifelse( x>x0, 1, ifelse( x<0, 0, outstress ) )
  outstress <- ifelse( outstress>1, 1, ifelse( outstress<0, 0, outstress ) )

  return( outstress )

}
