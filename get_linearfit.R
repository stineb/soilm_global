get_linearfit <- function( df, monthly=FALSE ){
  
  require(dplyr)
  
  if (monthly){
    ## add date and MOY to dataframe nice_agg
    df <- df %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
    df <- df %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

    ## aggregate nice_agg to monthly values
    df <- df %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_mean = mean( soilm_mean, na.rm=TRUE ) )    
  }

  ##------------------------------------------------------------------------
  ## Determine maximum LUE reduction (mean fLUE in lowest bin) for each site (-> df_flue0)
  ##------------------------------------------------------------------------
  ## Bin values and get mean fLUE for soil moisture < 0.25 for each site (:= flue0)
  intervals <- seq(0, 1, 0.25)
  df$ininterval <- NULL
  df <- df %>% mutate( ininterval = cut( soilm_mean , breaks = intervals ) ) %>% group_by( mysitename, ininterval )
  df_flue0 <- df %>%  dplyr::summarise( flue0=mean( fvar, na.rm=TRUE ) ) %>% 
                      complete( ininterval, fill = list( flue0 = NA ) ) %>% 
                      dplyr::filter( ininterval=="(0,0.25]" )

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
  df_flue0 <- df_flue0 %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

  ##------------------------------------------------------------------------
  ## Fit linear model
  ##------------------------------------------------------------------------
  linmod <- lm( flue0 ~ meanalpha, data=df_flue0 )

  return( list( linmod=linmod, df_flue0=df_flue0 ) )

}