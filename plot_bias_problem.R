library(dplyr)
library(tidyr)
library(tibble)

myboxplot <- function( ... ){
	boxplot( ..., staplewex=0, whisklty=1, outline=FALSE )
}

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" )

## IMPORTANT: REMOVE DUPLICATE ROWS (were introduced because one date can below to multiple drought instances (within their before-after window) )
df_dday_agg    <- df_dday_agg    %>% select( -dday, -inst ) %>% unique()
df_dday_8d_agg <- df_dday_8d_agg %>% select( -dday, -inst ) %>% unique()


# df_dday_8d_agg <- filter( df_dday_8d_agg, mysitename!="US-Var")

##------------------------------------------------
## bin data
##------------------------------------------------
## define fLUE bins
nbins <- 5
max <- 1.0
binwidth <- max/nbins
bins <- seq( from=0, to=max, by=binwidth )

## define soil moisture bins
sbins <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0)
nsbins <- length(sbins)
sbins_plot <- seq( from=0, to=1, by=1.0/(length(sbins)-1) )

df_dday_8d_agg <- df_dday_8d_agg %>% mutate( infbin = cut( fvar, breaks = bins ), insbin = cut( soilm_mean, breaks = sbins ) )
 
##------------------------------------------------
## normalise ratios 
## to 1 for bin fvar = 1.0-1.1
##------------------------------------------------
df_binmedians <- df_dday_8d_agg %>%   group_by( infbin ) %>% 
                                      summarise( median_modis   = median( bias_modis,   na.rm = TRUE ),
                                                 median_vpm     = median( bias_vpm,     na.rm = TRUE ),
                                                 median_bess_v1 = median( bias_bess_v1, na.rm = TRUE ),
                                                 median_pmodel  = median( bias_pmodel,  na.rm = TRUE )
                                                )

df_dday_8d_agg_norm <- df_dday_8d_agg %>%  mutate( gpp_modis   = gpp_modis     / df_binmedians$median_modis[nbins],
                                                   gpp_vpm     = gpp_vpm      / df_binmedians$median_vpm[nbins], 
                                                   gpp_bess_v1 = gpp_bess_v1 / df_binmedians$median_bess_v1[nbins], 
                                                   gpp_pmodel  = gpp_pmodel / df_binmedians$median_pmodel[nbins]  
                                                   ) %>%
                                          mutate( bias_pmodel_diff  = gpp_pmodel  - gpp_obs,
                                                  bias_modis_diff   = gpp_modis   - gpp_obs,
                                                  bias_vpm_diff     = gpp_vpm     - gpp_obs,
                                                  bias_bess_v1_diff = gpp_bess_v1 - gpp_obs
                                                 )

##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
library(tibble)
df_dday_8d_agg_norm <- df_dday_8d_agg_norm %>% as_tibble() %>%
                       mutate(  
                                bias_modis_diff_corr   = gpp_modis   * fvar - gpp_obs,
                                bias_vpm_diff_corr     = gpp_vpm     * fvar - gpp_obs, 
                                bias_bess_v1_diff_corr = gpp_bess_v1 * fvar - gpp_obs, 
                                bias_pmodel_diff_corr  = gpp_pmodel  * fvar - gpp_obs,

                                bias_modis_diff_corr_I   = gpp_modis   * flue_est_I - gpp_obs,
                                bias_vpm_diff_corr_I     = gpp_vpm     * flue_est_I - gpp_obs, 
                                bias_bess_v1_diff_corr_I = gpp_bess_v1 * flue_est_I - gpp_obs, 
                                bias_pmodel_diff_corr_I  = gpp_pmodel  * flue_est_I - gpp_obs,

                                bias_modis_diff_corr_II   = gpp_modis   * flue_est_II - gpp_obs,
                                bias_vpm_diff_corr_II     = gpp_vpm     * flue_est_II - gpp_obs, 
                                bias_bess_v1_diff_corr_II = gpp_bess_v1 * flue_est_II - gpp_obs, 
                                bias_pmodel_diff_corr_II  = gpp_pmodel  * flue_est_II - gpp_obs,

                                ## flue_est_IV is based on separate grasslands and others
                                bias_modis_diff_corr_IV   = gpp_modis   * flue_est_IV - gpp_obs,
                                bias_vpm_diff_corr_IV     = gpp_vpm     * flue_est_IV - gpp_obs, 
                                bias_bess_v1_diff_corr_IV = gpp_bess_v1 * flue_est_IV - gpp_obs, 
                                bias_pmodel_diff_corr_IV  = gpp_pmodel  * flue_est_IV - gpp_obs,

                                bias_modis_diff_corr_III   = gpp_modis   * flue_est_III - gpp_obs,
                                bias_vpm_diff_corr_III     = gpp_vpm     * flue_est_III - gpp_obs, 
                                bias_bess_v1_diff_corr_III = gpp_bess_v1 * flue_est_III - gpp_obs, 
                                bias_pmodel_diff_corr_III  = gpp_pmodel  * flue_est_III - gpp_obs,

                                ## flue_est_V is based on exponential stress function
                                bias_modis_diff_corr_V   = gpp_modis   * flue_est_V - gpp_obs,
                                bias_vpm_diff_corr_V     = gpp_vpm     * flue_est_V - gpp_obs, 
                                bias_bess_v1_diff_corr_V = gpp_bess_v1 * flue_est_V - gpp_obs, 
                                bias_pmodel_diff_corr_V  = gpp_pmodel  * flue_est_V - gpp_obs
                              )

df_dday_8d_agg <- df_dday_8d_agg %>% as_tibble() %>%
                     mutate(  
                              bias_modis_diff_corr   = gpp_modis   * fvar - gpp_obs,
                              bias_vpm_diff_corr     = gpp_vpm     * fvar - gpp_obs, 
                              bias_bess_v1_diff_corr = gpp_bess_v1 * fvar - gpp_obs, 
                              bias_pmodel_diff_corr  = gpp_pmodel  * fvar - gpp_obs,

                              bias_modis_diff_corr_I   = gpp_modis   * flue_est_I - gpp_obs,
                              bias_vpm_diff_corr_I     = gpp_vpm     * flue_est_I - gpp_obs, 
                              bias_bess_v1_diff_corr_I = gpp_bess_v1 * flue_est_I - gpp_obs, 
                              bias_pmodel_diff_corr_I  = gpp_pmodel  * flue_est_I - gpp_obs,

                              bias_modis_diff_corr_II   = gpp_modis   * flue_est_II - gpp_obs,
                              bias_vpm_diff_corr_II     = gpp_vpm     * flue_est_II - gpp_obs, 
                              bias_bess_v1_diff_corr_II = gpp_bess_v1 * flue_est_II - gpp_obs, 
                              bias_pmodel_diff_corr_II  = gpp_pmodel  * flue_est_II - gpp_obs,

                              ## flue_est_IV is based on separate grasslands and others
                              bias_modis_diff_corr_IV   = gpp_modis   * flue_est_IV - gpp_obs,
                              bias_vpm_diff_corr_IV     = gpp_vpm     * flue_est_IV - gpp_obs, 
                              bias_bess_v1_diff_corr_IV = gpp_bess_v1 * flue_est_IV - gpp_obs, 
                              bias_pmodel_diff_corr_IV  = gpp_pmodel  * flue_est_IV - gpp_obs,

                              bias_modis_diff_corr_III   = gpp_modis   * flue_est_III - gpp_obs,
                              bias_vpm_diff_corr_III     = gpp_vpm     * flue_est_III - gpp_obs, 
                              bias_bess_v1_diff_corr_III = gpp_bess_v1 * flue_est_III - gpp_obs, 
                              bias_pmodel_diff_corr_III  = gpp_pmodel  * flue_est_III - gpp_obs,

                              ## flue_est_V is based on exponential stress function
                              bias_modis_diff_corr_V   = gpp_modis   * flue_est_V - gpp_obs,
                              bias_vpm_diff_corr_V     = gpp_vpm     * flue_est_V - gpp_obs, 
                              bias_bess_v1_diff_corr_V = gpp_bess_v1 * flue_est_V - gpp_obs, 
                              bias_pmodel_diff_corr_V  = gpp_pmodel  * flue_est_V - gpp_obs
                            )

##------------------------------------------------
## Treat all different model's data as one - gather data
##------------------------------------------------
tmp  <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff) %>% 
                                gather( model, bias_diff, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff )


##------------------------------------------------
## Final plot: bias (the problem)
## showing:
## - uncorrected, normalised, pooled models
##------------------------------------------------
par(las=1)

## pooled
myboxplot( bias_diff ~ infbin, data = tmp, ylim=c(-5,6.5), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )


## p-model
myboxplot( bias_pmodel_diff ~ infbin, data = df_dday_8d_agg_norm, ylim=c(-5,6.5), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "P-model, normalised", line=0.5, adj = 0, font=2 )

## modis
myboxplot( bias_modis_diff ~ infbin, data = df_dday_8d_agg_norm, ylim=c(-5,6.5), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "MODIS, normalised", line=0.5, adj = 0, font=2 )

## VPM
myboxplot( bias_vpm_diff ~ infbin, data = df_dday_8d_agg_norm, ylim=c(-5,6.5), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "VPM, normalised", line=0.5, adj = 0, font=2 )

## BESS
myboxplot( bias_bess_diff ~ infbin, data = df_dday_8d_agg_norm, ylim=c(-5,6.5), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "BESS, normalised", line=0.5, adj = 0, font=2 )


