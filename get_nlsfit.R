get_nlsfit <- function( df, monthly=FALSE ){

  require(minpack.lm)
  require(dplyr)

  source("calc_flue_est_alpha.R")

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
  df <- df %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

  if (monthly){
    ## add date and MOY to dataframe nice_agg
    df <- df %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
    df <- df %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

    ## aggregate nice_agg to monthly values
    df <- df %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_mean = mean( soilm_mean, na.rm=TRUE ) )    
  }

  ##------------------------------------
  ## Estimate fLUE using non-linear least squares on quadratic function
  ##------------------------------------
  nlsfit <- try( 
                nlsLM( 
                      fvar ~ calc_flue_est_alpha( soilm_mean, meanalpha, apar, bpar, cpar, dpar ),
                      data=df,
                      start=list( apar=0.11, bpar=0.89, cpar=0.125, dpar=0.75 ),
                      lower=c( apar=0.001, bpar=0.5, cpar=0.001, dpar=0.3 ),
                      upper=c( apar=0.3, bpar=1.5, cpar=0.3, dpar=1.0 ),
                      algorithm="port"
                      ) 
                )

  return( nlsfit )

}