library(dplyr)
library(readr)

##------------------------------------------------
## Data filtered to fLUE-drought-affected periods
##------------------------------------------------
## Load aligned aggregated data; prepared by 'execute_reshape_align_nn_fluxnet2015.R', calling 'reshape_align_nn_fluxnet2015.R'
load( "data/data_aligned_agg.Rdata" )

## IMPORTANT: REMOVE DUPLICATE ROWS (were introduced because one date can below to multiple drought instances (within their before-after window) )
df_dday_agg    <- df_dday_agg    %>% select( -dday, -inst ) %>% unique()
df_dday_8d_agg <- df_dday_8d_agg %>% select( -dday, -inst ) %>% unique()

## Select variables for which data can be made publicly accessible
df_dday_agg    <- df_dday_agg    %>% select(  site_id = mysitename, date, soilm_splash, flue = fvar, beta_a = flue_est_I, 
																							beta_b = flue_est_IV, beta_c = flue_est_III, alpha, gpp_pmodel
																							)
df_dday_8d_agg <- df_dday_8d_agg %>% select(  site_id = mysitename, date, soilm_splash, flue = fvar, beta_a = flue_est_I, 
																							beta_b = flue_est_IV, beta_c = flue_est_III, alpha, gpp_pmodel
																							)

if (!dir.exists("./data_openaccess")) system("mkdir data_openaccess")

write_csv( df_dday_agg, path = "data_openaccess/gpp_alg_daily_fluxnet_stocker18natgeo.csv" )
write_csv( df_dday_8d_agg, path = "data_openaccess/gpp_alg_8daily_fluxnet_stocker18natgeo.csv" )

##------------------------------------------------
## Full data with reliable soil moisture information
##------------------------------------------------
# Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( "data/nice_all_agg_lue_obs_evi.Rdata" )      # loads 'nice_agg'
load( "data/nice_all_8d_agg_lue_obs_evi.Rdata" )   # loads 'nice_8d'

successcodes <- read.csv( "successcodes.csv" )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Use only sites where NN method worked (i.e. that had clear and identifiable soil moisture limitation)
nice_agg     <- nice_agg     %>% filter( mysitename %in% do.sites )
nice_8d_agg  <- nice_8d_agg  %>% filter( mysitename %in% do.sites )

## Select variables for which data can be made publicly accessible
nice_agg    <- nice_agg    %>% select( site_id = mysitename, date      , gpp_pmodel, aet_splash = aet_pmodel, pet_splash = pet_pmodel )
nice_8d_agg <- nice_8d_agg %>% select( site_id = mysitename, date_start, gpp_pmodel, aet_splash = aet_pmodel, pet_splash = pet_pmodel )

write_csv( nice_agg, path = "data_openaccess/gpp_daily_fluxnet_stocker18natgeo.csv" )
write_csv( nice_8d_agg, path = "data_openaccess/gpp_8daily_fluxnet_stocker18natgeo.csv" )
