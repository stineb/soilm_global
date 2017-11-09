compl_df_flue_est <- function( df, linearfit, x0_fix=0.9 ){

  ## this requires columns in 'df' called 'soilm_mean' and 'meanalpha'

  require(dplyr)

  source("stress_quad_1sided.R")

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  if (is.null(df$meanalpha)){
    load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
    df <- df %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )
  }

  ##------------------------------------
  ## add estimated fLUE values to data frame
  ##------------------------------------
  ## Estimate fLUE based on linear fit between fLUE0 and mean-alpha (fixed "tie points")
  df <- df %>% mutate( flue_est = stress_quad_1sided_alpha( soilm_mean, meanalpha, x0_fix, coef(linearfit2$linmod)[["(Intercept)"]], coef(linearfit2$linmod)[["meanalpha"]] ) )

  return( df )

}