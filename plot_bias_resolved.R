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


df_dday_8d_agg <- filter( df_dday_8d_agg, mysitename!="US-Var")

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
## Quick plots I:
##------------------------------------------------
## Bias for each model separately
## pmodel
par(las=1)
myboxplot( bias_pmodel_diff ~ infbin, data=df_dday_8d_agg, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( bias_modis_diff ~ infbin, data=df_dday_8d_agg, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff ~ infbin, data=df_dday_8d_agg, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff ~ infbin, data=df_dday_8d_agg, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS", line=0.5, adj = 0, font=2 )
 
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
## Quick plots II:
##------------------------------------------------
## Bias for each model separately after "normalisation"
## pmodel
myboxplot( bias_pmodel_diff ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, normalised", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( bias_modis_diff ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, normalised", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, normalised", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, normalised", line=0.5, adj = 0, font=2 )

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
## Quick plots III:
##------------------------------------------------
## Bias for each model separately after "normalisation"

## Corrected by fLUE

## pmodel
myboxplot( bias_pmodel_diff_corr ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, norm., corr. fLUE", line=0.5, adj = 0, font=2 )

# library(LSD)
# with( df_dday_8d_agg_norm, heatscatter( fvar, bias_pmodel_diff, ylim = c(-6,6), xlim = c(0,1.2) ) )
# abline( h=0, lty=3 )

# with( df_dday_8d_agg_norm, heatscatter( fvar, bias_pmodel_diff_corr, ylim = c(-6,6), xlim = c(0,1.2) ) )
# abline( h=0, lty=3 )

# with( df_dday_8d_agg_norm, heatscatter( fvar, bias_pmodel_diff_corr_III, ylim = c(-6,6), xlim = c(0,1.2) ) )
# abline( h=0, lty=3 )

## MODIS
myboxplot( bias_modis_diff_corr ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, norm., corr. fLUE", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff_corr ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, norm., corr. fLUE", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff_corr ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, norm., corr. fLUE", line=0.5, adj = 0, font=2 )


## Corrected by estimated fLUE (following approach II)

## pmodel
myboxplot( bias_pmodel_diff_corr_II ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, norm., corr. II", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( bias_modis_diff_corr_II ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, norm., corr. II", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff_corr_II ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, norm., corr. II", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff_corr_II ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, norm., corr. II", line=0.5, adj = 0, font=2 )


## Corrected by estimated fLUE (following approach IV)

## pmodel
myboxplot( bias_pmodel_diff_corr_IV ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, norm., corr. IV", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( bias_modis_diff_corr_IV ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, norm., corr. IV", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff_corr_IV ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, norm., corr. IV", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff_corr_IV ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, norm., corr. IV", line=0.5, adj = 0, font=2 )

## Corrected by estimated fLUE (following approach V)

## pmodel
myboxplot( bias_pmodel_diff_corr_V ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, norm., corr. V", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( bias_modis_diff_corr_V ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, norm., corr. V", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( bias_vpm_diff_corr_V ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, norm., corr. V", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bias_bess_v1_diff_corr_V ~ infbin, data=df_dday_8d_agg_norm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, norm., corr. V", line=0.5, adj = 0, font=2 )


##------------------------------------------------
## Treat data from different stress functions as one - gather data
##------------------------------------------------
pmodel <- df_dday_8d_agg_norm %>% select( mysitename, date, infbin, insbin, bias_pmodel_diff_corr_I, bias_pmodel_diff_corr_IV, bias_pmodel_diff_corr_III ) %>%
                                  gather( fct, pmodel, bias_pmodel_diff_corr_I, bias_pmodel_diff_corr_IV, bias_pmodel_diff_corr_III )

modis <- df_dday_8d_agg_norm %>% select( mysitename, date, infbin, insbin, bias_modis_diff_corr_I, bias_modis_diff_corr_IV, bias_modis_diff_corr_III ) %>%
                                 gather( fct, modis, bias_modis_diff_corr_I, bias_modis_diff_corr_IV, bias_modis_diff_corr_III )

vpm <- df_dday_8d_agg_norm %>% select( mysitename, date, infbin, insbin, bias_vpm_diff_corr_I, bias_vpm_diff_corr_IV, bias_vpm_diff_corr_III ) %>%
                               gather( fct, vpm, bias_vpm_diff_corr_I, bias_vpm_diff_corr_IV, bias_vpm_diff_corr_III )

bess_v1 <- df_dday_8d_agg_norm %>% select( mysitename, date, infbin, insbin, bias_bess_v1_diff_corr_I, bias_bess_v1_diff_corr_IV, bias_bess_v1_diff_corr_III ) %>%
                                   gather( fct, bess_v1, bias_bess_v1_diff_corr_I, bias_bess_v1_diff_corr_IV, bias_bess_v1_diff_corr_III )

## Corrected by estimated fLUE (following approach II)

## pmodel
myboxplot( pmodel ~ infbin, data=pmodel, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "P-model, norm., corr. pooled", line=0.5, adj = 0, font=2 )

## MODIS
myboxplot( modis ~ infbin, data=modis, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "MODIS, norm., corr. pooled", line=0.5, adj = 0, font=2 )

## vpm
myboxplot( vpm ~ infbin, data=vpm, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "VPM, norm., corr. pooled", line=0.5, adj = 0, font=2 )

## bess_v1
myboxplot( bess_v1 ~ infbin, data=bess_v1, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6) )
abline( h=0, lty=3 )
mtext( "BESS, norm., corr. pooled", line=0.5, adj = 0, font=2 )


##------------------------------------------------
## Treat all different model's data as one - gather data
##------------------------------------------------
tmp  <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff) %>% 
                                gather( model, bias_diff, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff )

tmp0 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr, bias_vpm_diff_corr, bias_bess_v1_diff_corr, bias_pmodel_diff_corr ) %>% 
                                gather( model, bias_diff_corr, bias_modis_diff_corr, bias_vpm_diff_corr, bias_bess_v1_diff_corr, bias_pmodel_diff_corr )

tmp1 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_I, bias_vpm_diff_corr_I, bias_bess_v1_diff_corr_I, bias_pmodel_diff_corr_I ) %>% 
                                gather( model, bias_diff_corr_I, bias_modis_diff_corr_I, bias_vpm_diff_corr_I, bias_bess_v1_diff_corr_I, bias_pmodel_diff_corr_I )

tmp2 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV ) %>% 
                                gather( model, bias_diff_corr_II, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV )

tmp3 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_III, bias_vpm_diff_corr_III, bias_bess_v1_diff_corr_III, bias_pmodel_diff_corr_III ) %>% 
                                gather( model, bias_diff_corr_III, bias_modis_diff_corr_III, bias_vpm_diff_corr_III, bias_bess_v1_diff_corr_III, bias_pmodel_diff_corr_III )

tmp4 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV ) %>% 
                                gather( model, bias_diff_corr_IV, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV )

tmp5 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_V, bias_vpm_diff_corr_V, bias_bess_v1_diff_corr_V, bias_pmodel_diff_corr_V ) %>% 
                                gather( model, bias_diff_corr_V, bias_modis_diff_corr_V, bias_vpm_diff_corr_V, bias_bess_v1_diff_corr_V, bias_pmodel_diff_corr_V )

##------------------------------------------------
## Quick plots IV:
##------------------------------------------------
par(las=1)
myboxplot( bias_diff ~ infbin, data = tmp, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr ~ infbin, data = tmp0, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. fLUE", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr_I ~ infbin, data = tmp1, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. I", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr_II ~ infbin, data = tmp2, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. II", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr_III ~ infbin, data = tmp3, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. III", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr_IV ~ infbin, data = tmp4, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. IV", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr_V ~ infbin, data = tmp5, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. V", line=0.5, adj = 0, font=2 )

pooled <- select( pmodel, -fct, bias=pmodel ) %>% 
          bind_rows( select( modis, -fct, bias=modis) ) %>% 
          bind_rows( select( vpm, -fct, bias=vpm) ) %>% 
          bind_rows( select( bess_v1, -fct, bias=bess_v1) )

##------------------------------------------------
## Quick plots V:
##------------------------------------------------
par(las=1)
myboxplot( bias_diff ~ infbin, data = tmp, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )

myboxplot( bias_diff_corr ~ infbin, data = tmp0, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., corr. fLUE", line=0.5, adj = 0, font=2 )

myboxplot( bias ~ infbin, data = pooled, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ) )
abline( h=0, lty=3 )
mtext( "Pooled models, norm., pooled corr.", line=0.5, adj = 0, font=2 )



##------------------------------------------------
## Final plot, pooled models, corrected by fLUE
## showing:
## - uncorrected, normalised, pooled models
## - corrected by I, IV, and III, normalised, pooled models
## - corrected by fLUE, normalised, pooled models
##------------------------------------------------
pdf("fig/bias_pooledmodels.pdf", width = 7, height = 6)
  par(las=1)
  space <- 0.1
  soff <- 0.017
  off <- 0.0475
              
  ## uncorrected
  myboxplot( bias_diff ~ infbin, data = tmp, ylim=c(-6,6), at=(bins[1:nbins]+1/nbins*space-0.02), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato", boxwex=0.05, xlim=c(-0.07,0.9) )

  ## corrected by fLUE
  myboxplot( bias_diff_corr ~ infbin, data = tmp0, add=TRUE, at=(bins[1:nbins]+1/nbins*space+0.03), col="royalblue3", boxwex=0.05, axes=FALSE )

  # ## corrected by IV
  # myboxplot( bias_diff_corr_IV ~ infbin, data = tmp4, add=TRUE, at=(bins[1:nbins]+1/nbins*space+soff*2+off), col="springgreen3", boxwex=0.05, axes=FALSE )

  abline( h=0, lty=3 )
  mtext( "", line=0.5, adj = 0, font=2 )
dev.off()
system("open fig/bias_pooledmodels.pdf")

##------------------------------------------------
## Final plot, P-model, corrected by fLUE and I, IV, III
## showing:
## - uncorrected, normalised, pooled models
## - corrected by I, IV, and III, normalised, pooled models
## - corrected by fLUE, normalised, pooled models
##------------------------------------------------
pdf("fig/bias_pmodel_resolved.pdf", width = 7, height = 6)
  par(las=1)
  space <- 0.1
  soff <- 0.017
  off <- 0.0475
              
  ## uncorrected
  par(las=1)
  myboxplot( bias_pmodel_diff ~ infbin, data=df_dday_8d_agg, at=(bins[1:nbins]+1/nbins*space-0.02), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), ylim=c(-6,6), col="tomato", boxwex=0.05, xlim=c(-0.05,0.95) )

  ## corrected by fLUE
  myboxplot( bias_pmodel_diff_corr ~ infbin, data=df_dday_8d_agg, at=(bins[1:nbins]+1/nbins*space+0.03), add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.05 )

  ## corrected by I
  myboxplot( bias_pmodel_diff_corr_I ~ infbin, data=df_dday_8d_agg, at=(bins[1:nbins]+1/nbins*space+soff*1+0.05), add=TRUE, col="springgreen1", boxwex=0.016, axes=FALSE )

  ## corrected by IV
  myboxplot( bias_pmodel_diff_corr_IV ~ infbin, data=df_dday_8d_agg, at=(bins[1:nbins]+1/nbins*space+soff*2+0.05), add=TRUE, col="springgreen3", boxwex=0.016, axes=FALSE )

  ## corrected by III
  myboxplot( bias_pmodel_diff_corr_III ~ infbin, data=df_dday_8d_agg, at=(bins[1:nbins]+1/nbins*space+soff*3+0.05), add=TRUE, col="springgreen4", boxwex=0.016, axes=FALSE )

  abline( h=0, lty=3 )
  mtext( "", line=0.5, adj = 0, font=2 )
dev.off()
system("open fig/bias_pmodel_resolved.pdf")


##------------------------------------------------
## Final plot: bias (the problem)
## showing:
## - uncorrected, normalised, pooled models
##------------------------------------------------
par(las=1)
myboxplot( bias_diff ~ infbin, data = tmp, ylim=c(-7,7), xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), col="tomato" )
abline( h=0, lty=3 )
mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )

