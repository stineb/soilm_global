get_linearfit3 <- function( df, useweights=FALSE ){
  ##------------------------------------------------------------------------
  ## Fit a and b parameters to minimise the difference between GPPobs and
  ## GPP_Pmodel * fLUEest
  ## This is different from linearfit2 because it fits directly.
  ##------------------------------------------------------------------------
  require(minpack.lm)

  source("stress_quad_1sided.R")

  if (useweights){
  	weights <- 1.0 - df$fvar 
    weights <- ifelse( weights<0, 0, weights )
    weights <- ifelse( is.na(weights), 0, weights )
    # weights <- weights^2
  } else {
    weights <- rep( 1.0, nrow(df) )
  }

	fit <- try(
							nlsLM(
										ratio_obs_mod_pmodel ~ stress_quad_1sided_alpha( soilm_splash220, meanalpha, x0=0.9, apar, bpar ),
										data = df,
										start = list( apar=0.2, bpar=0.5 ),
										lower = c( apar=-1, bpar=0.0 ),
										upper = c( apar=0.5, bpar=2.0 ),
                    algorithm="port", 
                    weights = weights
										) 
							)
	return( fit )
}