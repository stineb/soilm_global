get_nlsfit <- function( df, target="fvar", monthly=FALSE, bin=FALSE, x0_fix=0.9 ){

  require(minpack.lm)
  require(dplyr)
  require(tidyr)

  source("calc_flue_est_alpha.R")
  source("stress_quad_1sided.R")

  if (monthly){
    ## add date and MOY to dataframe nice_agg
    df <- df %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
    df <- df %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

    ## aggregate nice_agg to monthly values
    df <- df %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_mean = mean( soilm_mean, na.rm=TRUE ) )    
  }

  if (bin){
    nsbins <- 10
    sbins <- seq( 0.0, 1.0, 1.0/nsbins )
    xvals <- sbins[1:(length(sbins)-1)] + (sbins[2]-sbins[1])/2

    ## Get median by interval and get fvar_vs_soilm for this site (used for clustering)
    df$insbin <- NULL
    df <- df %>% mutate( insbin = cut( soilm_mean, breaks = sbins ) ) %>% group_by( insbin, mysitename )
    df <- df %>%  summarise( soilm_mean=mean( soilm_mean, na.rm=TRUE ), fvar=mean( fvar, na.rm=TRUE ), ratio_obs_mod=mean( ratio_obs_mod, na.rm=TRUE ) ) %>% 
                  complete( insbin, fill = list( fvar = NA ) )
  }

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  if (is.null(df$meanalpha)){
    load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
    df <- df %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )
  }
  
  ##------------------------------------
  ## Estimate fLUE using non-linear least squares on quadratic function
  ##------------------------------------
  eq <- paste0( target, " ~ stress_quad_1sided_alpha( soilm_mean, meanalpha, x0, apar, bpar )" )
  nlsfit <- try( 
                nlsLM( 
                      # fvar ~ calc_flue_est_alpha( soilm_mean, meanalpha, apar, bpar, cpar, dpar ),
                      # fvar ~ stress_quad_1sided_alpha( soilm_mean, meanalpha, x0, apar, bpar ),
                      eq,
                      data=df,
                      start=list( x0=0.9, apar=-0.19,  bpar=0.96 ),
                      lower=c(       0.9, apar=-1   ,  bpar=0.1  ),
                      upper=c(       0.9, apar=1    ,  bpar=2    ),
                      algorithm="port"
                      ) 
                )

  return( nlsfit )

}