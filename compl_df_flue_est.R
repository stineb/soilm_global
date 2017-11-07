compl_df_flue_est <- function( df, linearfit, nlsfit ){

  ## this requires columns in 'df' called 'soilm_mean' and 'meanalpha'

  require(dplyr)

  source("calc_flue_est_alpha.R")

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  if (is.null(df$meanalpha)){
    load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
    df <- df %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )
  }

  ##------------------------------------
  ## add estimated fLUE values to data frame
  ##------------------------------------
  ## Estimate fLUE based on linear fit between fLUE0 and mean-alpha (fixed "tie points")
  df <- df %>% mutate( flue_est = calc_flue_est_alpha( soilm_mean, meanalpha, coef(linearfit$linmod)[1], coef(linearfit$linmod)[2], 0.125, 0.75 ) )

  ## Estimate fLUE based on fully fitted relationship using non-linear least squares fit function
  df <- df %>% mutate( flue_est_nls = calc_flue_est_alpha( soilm_mean, meanalpha, coef(nlsfit)[[ "apar" ]], coef(nlsfit)[[ "bpar" ]], coef(nlsfit)[[ "cpar" ]], coef(nlsfit)[[ "dpar" ]] ) )

  return( df )

}