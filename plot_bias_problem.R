library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)

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
sbins <- c(0, 0.2, 0.3, 0.4, 0.6, 1.0)
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
                                                 median_pmodel  = median( bias_pmodel,  na.rm = TRUE ),
                                                 median_mte     = median( bias_mte   ,  na.rm = TRUE )
                                                )

df_dday_8d_agg_norm <- df_dday_8d_agg %>%  mutate( gpp_modis   = gpp_modis     / df_binmedians$median_modis[nbins],
                                                   gpp_vpm     = gpp_vpm      / df_binmedians$median_vpm[nbins], 
                                                   gpp_bess_v1 = gpp_bess_v1 / df_binmedians$median_bess_v1[nbins], 
                                                   gpp_pmodel  = gpp_pmodel / df_binmedians$median_pmodel[nbins],
                                                   gpp_mte     = gpp_mte   / df_binmedians$median_mte[nbins]  
                                                   ) %>%
                                          mutate( bias_pmodel_diff  = gpp_pmodel  - gpp_obs,
                                                  bias_modis_diff   = gpp_modis   - gpp_obs,
                                                  bias_vpm_diff     = gpp_vpm     - gpp_obs,
                                                  bias_bess_v1_diff = gpp_bess_v1 - gpp_obs,
                                                  bias_mte_diff     = gpp_mte     - gpp_obs
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
                                bias_mte_diff_corr     = gpp_mte     * fvar - gpp_obs,

                                bias_modis_diff_corr_I   = gpp_modis   * flue_est_I - gpp_obs,
                                bias_vpm_diff_corr_I     = gpp_vpm     * flue_est_I - gpp_obs, 
                                bias_bess_v1_diff_corr_I = gpp_bess_v1 * flue_est_I - gpp_obs, 
                                bias_pmodel_diff_corr_I  = gpp_pmodel  * flue_est_I - gpp_obs,
                                bias_mte_diff_corr_I     = gpp_mte     * flue_est_I - gpp_obs,

                                bias_modis_diff_corr_II   = gpp_modis   * flue_est_II - gpp_obs,
                                bias_vpm_diff_corr_II     = gpp_vpm     * flue_est_II - gpp_obs, 
                                bias_bess_v1_diff_corr_II = gpp_bess_v1 * flue_est_II - gpp_obs, 
                                bias_pmodel_diff_corr_II  = gpp_pmodel  * flue_est_II - gpp_obs,
                                bias_mte_diff_corr_II     = gpp_mte     * flue_est_II - gpp_obs,

                                ## flue_est_IV is based on separate grasslands and others
                                bias_modis_diff_corr_IV   = gpp_modis   * flue_est_IV - gpp_obs,
                                bias_vpm_diff_corr_IV     = gpp_vpm     * flue_est_IV - gpp_obs, 
                                bias_bess_v1_diff_corr_IV = gpp_bess_v1 * flue_est_IV - gpp_obs, 
                                bias_pmodel_diff_corr_IV  = gpp_pmodel  * flue_est_IV - gpp_obs,
                                bias_mte_diff_corr_IV     = gpp_mte     * flue_est_IV - gpp_obs,

                                bias_modis_diff_corr_III   = gpp_modis   * flue_est_III - gpp_obs,
                                bias_vpm_diff_corr_III     = gpp_vpm     * flue_est_III - gpp_obs, 
                                bias_bess_v1_diff_corr_III = gpp_bess_v1 * flue_est_III - gpp_obs, 
                                bias_pmodel_diff_corr_III  = gpp_pmodel  * flue_est_III - gpp_obs,
                                bias_mte_diff_corr_III     = gpp_mte     * flue_est_III - gpp_obs,

                                ## flue_est_V is based on exponential stress function
                                bias_modis_diff_corr_V   = gpp_modis   * flue_est_V - gpp_obs,
                                bias_vpm_diff_corr_V     = gpp_vpm     * flue_est_V - gpp_obs, 
                                bias_bess_v1_diff_corr_V = gpp_bess_v1 * flue_est_V - gpp_obs, 
                                bias_pmodel_diff_corr_V  = gpp_pmodel  * flue_est_V - gpp_obs,
                                bias_mte_diff_corr_V     = gpp_mte     * flue_est_V - gpp_obs
                              )

df_dday_8d_agg <- df_dday_8d_agg %>% as_tibble() %>%
                     mutate(  
                              bias_modis_diff_corr   = gpp_modis   * fvar - gpp_obs,
                              bias_vpm_diff_corr     = gpp_vpm     * fvar - gpp_obs, 
                              bias_bess_v1_diff_corr = gpp_bess_v1 * fvar - gpp_obs, 
                              bias_pmodel_diff_corr  = gpp_pmodel  * fvar - gpp_obs,
                              bias_mte_diff_corr     = gpp_mte     * fvar - gpp_obs,

                              bias_modis_diff_corr_I   = gpp_modis   * flue_est_I - gpp_obs,
                              bias_vpm_diff_corr_I     = gpp_vpm     * flue_est_I - gpp_obs, 
                              bias_bess_v1_diff_corr_I = gpp_bess_v1 * flue_est_I - gpp_obs, 
                              bias_pmodel_diff_corr_I  = gpp_pmodel  * flue_est_I - gpp_obs,
                              bias_mte_diff_corr_I     = gpp_mte     * flue_est_I - gpp_obs,

                              bias_modis_diff_corr_II   = gpp_modis   * flue_est_II - gpp_obs,
                              bias_vpm_diff_corr_II     = gpp_vpm     * flue_est_II - gpp_obs, 
                              bias_bess_v1_diff_corr_II = gpp_bess_v1 * flue_est_II - gpp_obs, 
                              bias_pmodel_diff_corr_II  = gpp_pmodel  * flue_est_II - gpp_obs,
                              bias_mte_diff_corr_II     = gpp_mte     * flue_est_II - gpp_obs,

                              ## flue_est_IV is based on separate grasslands and others
                              bias_modis_diff_corr_IV   = gpp_modis   * flue_est_IV - gpp_obs,
                              bias_vpm_diff_corr_IV     = gpp_vpm     * flue_est_IV - gpp_obs, 
                              bias_bess_v1_diff_corr_IV = gpp_bess_v1 * flue_est_IV - gpp_obs, 
                              bias_pmodel_diff_corr_IV  = gpp_pmodel  * flue_est_IV - gpp_obs,
                              bias_mte_diff_corr_IV     = gpp_mte     * flue_est_IV - gpp_obs,

                              bias_modis_diff_corr_III   = gpp_modis   * flue_est_III - gpp_obs,
                              bias_vpm_diff_corr_III     = gpp_vpm     * flue_est_III - gpp_obs, 
                              bias_bess_v1_diff_corr_III = gpp_bess_v1 * flue_est_III - gpp_obs, 
                              bias_pmodel_diff_corr_III  = gpp_pmodel  * flue_est_III - gpp_obs,
                              bias_mte_diff_corr_III     = gpp_mte     * flue_est_III - gpp_obs,

                              ## flue_est_V is based on exponential stress function
                              bias_modis_diff_corr_V   = gpp_modis   * flue_est_V - gpp_obs,
                              bias_vpm_diff_corr_V     = gpp_vpm     * flue_est_V - gpp_obs, 
                              bias_bess_v1_diff_corr_V = gpp_bess_v1 * flue_est_V - gpp_obs, 
                              bias_pmodel_diff_corr_V  = gpp_pmodel  * flue_est_V - gpp_obs,
                              bias_mte_diff_corr_V     = gpp_mte     * flue_est_V - gpp_obs
                            )

##------------------------------------------------
## Treat all different model's data as one - gather data
##------------------------------------------------
tmp  <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff, bias_mte_diff) %>% 
                                gather( model, bias_diff, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff ) # , bias_mte_diff


##------------------------------------------------
## Final plot: bias (the problem)
## showing bias for:
## - uncorrected, normalised, pooled models
##------------------------------------------------
xlim <- c(0.5,5.5)
ylim <- c(-5,6.5)

## pooled, fLUE bins
pdf("fig/bias_pooled.pdf", width = 7, height = 6)
  par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
  plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
  rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
  myboxplot( bias_diff ~ infbin, data = tmp, col="tomato", add=TRUE )
  abline( h=0, lty=3 )
  mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )
dev.off()

## pooled, soil moisture bins
pdf("fig/bias_pooled_smbin.pdf", width = 7, height = 6)
  par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
  plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "soil moisture bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
  rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
  myboxplot( bias_diff ~ insbin, data = tmp, col="tomato", add=TRUE )
  abline( h=0, lty=3 )
  mtext( "Pooled models, normalised", line=0.5, adj = 0, font=2 )
dev.off()

## by individual models
pdf("fig/bias_bymodels_norm.pdf", width = 7, height = 6)
  par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
  plot( xlim, ylim, type="n", ylim=ylim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), xlim=xlim, axes=FALSE )
  rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
  myboxplot( bias_pmodel_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.3, add=TRUE, col="tomato", boxwex=0.2 )

  ## modis
  myboxplot( bias_modis_diff ~ infbin, data = df_dday_8d_agg_norm, add=TRUE, at=1:5-0.1, col="orchid", axes=FALSE, boxwex=0.2 )

  ## VPM
  myboxplot( bias_vpm_diff ~ infbin, data = df_dday_8d_agg_norm, add=TRUE, at=1:5+0.1, col="darkgoldenrod1", axes=FALSE, boxwex=0.2 )

  ## BESS
  myboxplot( bias_bess_v1_diff ~ infbin, data = df_dday_8d_agg_norm, add=TRUE, at=1:5+0.3, col="springgreen", axes=FALSE, boxwex=0.2 )

  # ## MTE
  # myboxplot( bias_mte_diff ~ infbin, data = df_dday_8d_agg_norm, add=TRUE, at=1:5-0.4, col="tomato", axes=FALSE, boxwex=0.2 )

  abline( h=0, lty=3 )
  legend("bottomleft", c("P-model", "MODIS", "VPM", "BESS"), fill=c("tomato", "orchid", "darkgoldenrod1", "springgreen"), bty="n")
dev.off()


## NOT NORMALISED, by individual models
pdf("fig/bias_bymodels.pdf", width = 7, height = 6)
  par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
  plot( xlim, ylim, type="n", ylim=ylim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), xlim=xlim, axes=FALSE )
  rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
  myboxplot( bias_pmodel_diff ~ infbin, data = df_dday_8d_agg, at=1:5-0.3, add=TRUE, col="tomato", boxwex=0.2 )

  ## modis
  myboxplot( bias_modis_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5-0.1, col="orchid", axes=FALSE, boxwex=0.2 )

  ## VPM
  myboxplot( bias_vpm_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5+0.1, col="darkgoldenrod1", axes=FALSE, boxwex=0.2 )

  ## BESS
  myboxplot( bias_bess_v1_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5+0.3, col="springgreen", axes=FALSE, boxwex=0.2 )

  # ## MTE
  # myboxplot( bias_mte_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5-0.4, col="tomato", axes=FALSE, boxwex=0.2 )

  abline( h=0, lty=3 )
  legend("bottomleft", c("P-model", "MODIS", "VPM", "BESS"), fill=c("tomato", "orchid", "darkgoldenrod1", "springgreen"), bty="n")
dev.off()


##------------------------------------------------
## Final plot: bias (the problem)
## showing ratio obs/mod for:
##------------------------------------------------
xlim <- c(0.5,5.5)
ylim <- c(0,4)

## by individual models
pdf("fig/ratio_obs_mod_bymodels.pdf", width = 7, height = 6)
  par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
  plot( xlim, ylim, type="n", ylim=ylim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), xlim=xlim, axes=FALSE )
  rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
  myboxplot( ratio_obs_mod_pmodel ~ infbin, data = df_dday_8d_agg, at=1:5-0.3, add=TRUE, col="tomato", boxwex=0.2 )

  ## modis
  myboxplot( ratio_obs_mod_modis ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5-0.1, col="orchid", axes=FALSE, boxwex=0.2 )

  ## VPM
  myboxplot( ratio_obs_mod_vpm ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5+0.1, col="darkgoldenrod1", axes=FALSE, boxwex=0.2 )

  ## BESS
  myboxplot( ratio_obs_mod_bess_v1 ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5+0.3, col="springgreen", axes=FALSE, boxwex=0.2 )

  # ## MTE
  # myboxplot( bias_mte_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=1:5-0.4, col="tomato", axes=FALSE, boxwex=0.2 )

  abline( h=1, lty=3 )
  legend("topleft", c("P-model", "MODIS", "VPM", "BESS"), fill=c("tomato", "orchid", "darkgoldenrod1", "springgreen"), bty="n")
dev.off()

