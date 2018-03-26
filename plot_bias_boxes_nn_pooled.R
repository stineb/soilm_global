require(dplyr)
require(tidyr)

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_8d_agg', 'df_dday_mte_agg', 'df_dday_bess_agg', 'df_dday_vpm_agg'

##------------------------------------------------
## bin data
##------------------------------------------------
## define fLUE bins
nbins <- 11
max <- 1.1
binwidth <- max/nbins
fvarbins <- seq( from=0, to=max, by=binwidth )
xvals <- fvarbins[1:nbins]+binwidth/2

df_dday_8d_agg <- df_dday_8d_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )

##------------------------------------------------
## Quick plots I:
##------------------------------------------------
## Bias for each model separately
## pmodel
par(las=1)
boxplot( ratio_obs_mod_pmodel ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "P-model", line=0.5, adj = 0, font=2 )

## MODIS
boxplot( ratio_obs_mod_modis ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "MODIS", line=0.5, adj = 0, font=2 )

## vpm
boxplot( ratio_obs_mod_vpm ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "VPM", line=0.5, adj = 0, font=2 )

## bess_v1
boxplot( ratio_obs_mod_bess_v1 ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "BESS", line=0.5, adj = 0, font=2 )

##------------------------------------------------
## normalise ratios 
## to 1 for bin fvar = 1.0-1.1
##------------------------------------------------
df_binmedians <- df_dday_8d_agg %>%   group_by( inbin ) %>% 
                                      summarise( median_modis   = median( ratio_obs_mod_modis,   na.rm = TRUE ),
                                                 median_vpm     = median( ratio_obs_mod_vpm,     na.rm = TRUE ),
                                                 median_bess_v1 = median( ratio_obs_mod_bess_v1, na.rm = TRUE ),
                                                 median_pmodel  = median( ratio_obs_mod_pmodel,  na.rm = TRUE )
                                                )

df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  ratio_obs_mod_modis   = ratio_obs_mod_modis     / df_binmedians$median_modis[nbins],
                                              ratio_obs_mod_vpm     = ratio_obs_mod_vpm      / df_binmedians$median_vpm[nbins], 
                                              ratio_obs_mod_bess_v1 = ratio_obs_mod_bess_v1 / df_binmedians$median_bess_v1[nbins], 
                                              ratio_obs_mod_pmodel  = ratio_obs_mod_pmodel / df_binmedians$median_pmodel[nbins]  
                                              )

##------------------------------------------------
## Quick plots II:
##------------------------------------------------
## Bias for each model separately after "normalisation"
## pmodel
boxplot( ratio_obs_mod_pmodel ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "P-model, normalised", line=0.5, adj = 0, font=2 )

## MODIS
boxplot( ratio_obs_mod_modis ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "MODIS, normalised", line=0.5, adj = 0, font=2 )

## vpm
boxplot( ratio_obs_mod_vpm ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "VPM, normalised", line=0.5, adj = 0, font=2 )

## bess_v1
boxplot( ratio_obs_mod_bess_v1 ~ inbin, data=filter(df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "BESS, normalised", line=0.5, adj = 0, font=2 )


##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
df_dday_8d_agg <- df_dday_8d_agg %>% as_tibble() %>%
                                     mutate(  
                                              ratio_obs_mod_modis_corr_1   = ifelse( flue_est_1>0, ratio_obs_mod_modis     / flue_est_1, NA ),
                                              ratio_obs_mod_vpm_corr_1     = ifelse( flue_est_1>0, ratio_obs_mod_vpm      / flue_est_1, NA ), 
                                              ratio_obs_mod_bess_v1_corr_1 = ifelse( flue_est_1>0, ratio_obs_mod_bess_v1 / flue_est_1, NA ), 
                                              ratio_obs_mod_pmodel_corr_1  = ifelse( flue_est_1>0, ratio_obs_mod_pmodel / flue_est_1, NA ),

                                              ratio_obs_mod_modis_corr_2   = ifelse( flue_est_2>0, ratio_obs_mod_modis     / flue_est_2, NA ),
                                              ratio_obs_mod_vpm_corr_2     = ifelse( flue_est_2>0, ratio_obs_mod_vpm      / flue_est_2, NA ), 
                                              ratio_obs_mod_bess_v1_corr_2 = ifelse( flue_est_2>0, ratio_obs_mod_bess_v1 / flue_est_2, NA ), 
                                              ratio_obs_mod_pmodel_corr_2  = ifelse( flue_est_2>0, ratio_obs_mod_pmodel / flue_est_2, NA ),

                                              ratio_obs_mod_modis_corr_3   = ifelse( flue_est_3>0, ratio_obs_mod_modis     / flue_est_3, NA ),
                                              ratio_obs_mod_vpm_corr_3     = ifelse( flue_est_3>0, ratio_obs_mod_vpm      / flue_est_3, NA ), 
                                              ratio_obs_mod_bess_v1_corr_3 = ifelse( flue_est_3>0, ratio_obs_mod_bess_v1 / flue_est_3, NA ), 
                                              ratio_obs_mod_pmodel_corr_3  = ifelse( flue_est_3>0, ratio_obs_mod_pmodel / flue_est_3, NA )
                                            ) %>%

                                     select( mysitename, fvar, date, dday, inst, inbin,
                                             ratio_obs_mod_modis, ratio_obs_mod_vpm, ratio_obs_mod_bess_v1, ratio_obs_mod_pmodel, 
                                             ratio_obs_mod_modis_corr_1, ratio_obs_mod_vpm_corr_1, ratio_obs_mod_bess_v1_corr_1, ratio_obs_mod_pmodel_corr_1, 
                                             ratio_obs_mod_modis_corr_2, ratio_obs_mod_vpm_corr_2, ratio_obs_mod_bess_v1_corr_2, ratio_obs_mod_pmodel_corr_2, 
                                             ratio_obs_mod_modis_corr_3, ratio_obs_mod_vpm_corr_3, ratio_obs_mod_bess_v1_corr_3, ratio_obs_mod_pmodel_corr_3 
                                             )

##------------------------------------------------
## Quick plots III:
##------------------------------------------------
## Bias for each model separately after "normalisation"

## approach II

## pmodel
boxplot( ratio_obs_mod_pmodel_corr_2 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "P-model, norm., corr. II", line=0.5, adj = 0, font=2 )

## MODIS
boxplot( ratio_obs_mod_modis_corr_2 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "MODIS, norm., corr. II", line=0.5, adj = 0, font=2 )

## vpm
boxplot( ratio_obs_mod_vpm_corr_2 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "VPM, norm., corr. II", line=0.5, adj = 0, font=2 )

## bess_v1
boxplot( ratio_obs_mod_bess_v1_corr_2 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "BESS, norm., corr. II", line=0.5, adj = 0, font=2 )


## approach III

## pmodel
boxplot( ratio_obs_mod_pmodel_corr_3 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "P-model, norm., corr. III", line=0.5, adj = 0, font=2 )

## MODIS
boxplot( ratio_obs_mod_modis_corr_3 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "MODIS, norm., corr. III", line=0.5, adj = 0, font=2 )

## vpm
boxplot( ratio_obs_mod_vpm_corr_3 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "VPM, norm., corr. III", line=0.5, adj = 0, font=2 )

## bess_v1
boxplot( ratio_obs_mod_bess_v1_corr_3 ~ inbin, data=filter( df_dday_8d_agg, mysitename!="US-Var"), xlab = "fLUE bin", ylab = "ratio obs./mod.", outline=FALSE, ylim=c(0,3) )
abline( h=1, lty=3 )
mtext( "BESS, norm., corr. III", line=0.5, adj = 0, font=2 )


##------------------------------------------------
## Treat all different model's data as one - gather data
##------------------------------------------------
tmp0 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_modis, ratio_obs_mod_vpm, ratio_obs_mod_bess_v1, ratio_obs_mod_pmodel) %>% 
                            gather( model, ratio_obs_mod, ratio_obs_mod_modis, ratio_obs_mod_vpm, ratio_obs_mod_bess_v1, ratio_obs_mod_pmodel )

tmp1 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_modis_corr_1, ratio_obs_mod_vpm_corr_1, ratio_obs_mod_bess_v1_corr_1, ratio_obs_mod_pmodel_corr_1 ) %>% 
                            gather( model, ratio_obs_mod_corr_1, ratio_obs_mod_modis_corr_1, ratio_obs_mod_vpm_corr_1, ratio_obs_mod_bess_v1_corr_1, ratio_obs_mod_pmodel_corr_1 )

tmp2 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_modis_corr_2, ratio_obs_mod_vpm_corr_2, ratio_obs_mod_bess_v1_corr_2, ratio_obs_mod_pmodel_corr_2 ) %>% 
                            gather( model, ratio_obs_mod_corr_2, ratio_obs_mod_modis_corr_2, ratio_obs_mod_vpm_corr_2, ratio_obs_mod_bess_v1_corr_2, ratio_obs_mod_pmodel_corr_2 )

tmp3 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_modis_corr_3, ratio_obs_mod_vpm_corr_3, ratio_obs_mod_bess_v1_corr_3, ratio_obs_mod_pmodel_corr_3 ) %>% 
                            gather( model, ratio_obs_mod_corr_3, ratio_obs_mod_modis_corr_3, ratio_obs_mod_vpm_corr_3, ratio_obs_mod_bess_v1_corr_3, ratio_obs_mod_pmodel_corr_3 )

df_dday_8d_agg_pooled <- tmp0 %>% left_join( tmp1, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") ) %>% 
                                  left_join( tmp2, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") ) %>% 
                                  left_join( tmp3, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") )

##------------------------------------------------
## Quick plots IV:
##------------------------------------------------
boxplot( ratio_obs_mod ~ inbin, data = filter( df_dday_8d_agg_pooled, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_corr_1 ~ inbin, data = filter( df_dday_8d_agg_pooled, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, norm., corr. I", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_corr_2 ~ inbin, data = filter( df_dday_8d_agg_pooled, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, norm., corr. II", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_corr_3 ~ inbin, data = filter( df_dday_8d_agg_pooled, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, norm., corr. III", line=0.5, adj = 0, font=2 )

##------------------------------------------------
## Treat all correction's data as one - gather data
##------------------------------------------------
tmp1 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_modis_corr_1, ratio_obs_mod_modis_corr_2, ratio_obs_mod_modis_corr_3 ) %>% 
                            gather( approach, ratio_obs_mod_modis_corr, ratio_obs_mod_modis_corr_1, ratio_obs_mod_modis_corr_2, ratio_obs_mod_modis_corr_3 )

tmp2 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_vpm_corr_1, ratio_obs_mod_vpm_corr_2, ratio_obs_mod_vpm_corr_3 ) %>% 
                            gather( approach, ratio_obs_mod_vpm_corr, ratio_obs_mod_vpm_corr_1, ratio_obs_mod_vpm_corr_2, ratio_obs_mod_vpm_corr_3 )

tmp3 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_pmodel_corr_1, ratio_obs_mod_pmodel_corr_2, ratio_obs_mod_pmodel_corr_3 ) %>% 
                            gather( approach, ratio_obs_mod_pmodel_corr, ratio_obs_mod_pmodel_corr_1, ratio_obs_mod_pmodel_corr_2, ratio_obs_mod_pmodel_corr_3 )

tmp4 <- df_dday_8d_agg %>%  select( mysitename, date, dday, inst, fvar, inbin, ratio_obs_mod_bess_v1_corr_1, ratio_obs_mod_bess_v1_corr_2, ratio_obs_mod_bess_v1_corr_3 ) %>% 
                            gather( approach, ratio_obs_mod_bess_v1_corr, ratio_obs_mod_bess_v1_corr_1, ratio_obs_mod_bess_v1_corr_2, ratio_obs_mod_bess_v1_corr_3 )

df_dday_8d_agg_pooledapproach <- tmp1 %>% left_join( tmp2, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") ) %>% 
                                          left_join( tmp3, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") ) %>% 
                                          left_join( tmp4, by = c("mysitename", "date", "dday", "inst", "fvar", "inbin") )

##------------------------------------------------
## Quick plots V:
##------------------------------------------------
boxplot( ratio_obs_mod_pmodel_corr ~ inbin, data = filter( df_dday_8d_agg_pooledapproach, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "P-model, normalised, pooled corr.", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_modis_corr ~ inbin, data = filter( df_dday_8d_agg_pooledapproach, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "MODIS, normalised, pooled corr.", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_bess_v1_corr ~ inbin, data = filter( df_dday_8d_agg_pooledapproach, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "BESS, normalised, pooled corr.", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_vpm_corr ~ inbin, data = filter( df_dday_8d_agg_pooledapproach, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "VPM, normalised, pooled corr.", line=0.5, adj = 0, font=2 )


##------------------------------------------------
## Gather by correction function: ratios from all three soil moisture stress factors
##------------------------------------------------
df_dday_8d_agg_pooled2 <- df_dday_8d_agg_pooled %>% select( mysitename, inbin, ratio_obs_mod_corr_1, ratio_obs_mod_corr_2, ratio_obs_mod_corr_3 ) %>% 
                                                    gather( level, ratio_obs_mod_corr, ratio_obs_mod_corr_1, ratio_obs_mod_corr_2, ratio_obs_mod_corr_3 )

##------------------------------------------------
## Quick plots VI:
##------------------------------------------------
boxplot( ratio_obs_mod ~ inbin, data = filter( df_dday_8d_agg_pooled, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, norm.", line=0.5, adj = 0, font=2 )

boxplot( ratio_obs_mod_corr ~ inbin, data = filter( df_dday_8d_agg_pooled2, mysitename != "US-Var"), outline=FALSE, ylim=c(0,5), xlab = "fLUE bin", ylab = "ratio obs./mod." )
abline( h=1, lty=3 )
mtext( "Pooled models, norm., pooled corr.", line=0.5, adj = 0, font=2 )




