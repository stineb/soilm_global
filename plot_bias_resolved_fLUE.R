library(dplyr, quietly = TRUE, warn.conflicts = FALSE )
library(tidyr, quietly = TRUE, warn.conflicts = FALSE )
library(tibble, quietly = TRUE, warn.conflicts = FALSE )
library(hydroGOF, quietly = TRUE, warn.conflicts = FALSE )
library(tibble, quietly = TRUE, warn.conflicts = FALSE )

myboxplot <- function( ... ){
	boxplot( ..., staplewex=0, whisklty=1, outline=FALSE )
}

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" )
# df_dday_8d_agg <- filter( df_dday_8d_agg, mysitename!="US-Var")

get_bias_resolved_flue <- function( df_dday_8d_agg ){

  ## IMPORTANT: REMOVE DUPLICATE ROWS (were introduced because one date can below to multiple drought instances (within their before-after window) )
  df_dday_8d_agg <- df_dday_8d_agg %>% select( -dday, -inst ) %>% unique()

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
  tmp  <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff, gpp_obs ) %>% 
                                  gather( model, bias_diff, bias_modis_diff, bias_vpm_diff, bias_bess_v1_diff, bias_pmodel_diff )

  tmp0 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr, bias_vpm_diff_corr, bias_bess_v1_diff_corr, bias_pmodel_diff_corr, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr, bias_modis_diff_corr, bias_vpm_diff_corr, bias_bess_v1_diff_corr, bias_pmodel_diff_corr )

  tmp1 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_I, bias_vpm_diff_corr_I, bias_bess_v1_diff_corr_I, bias_pmodel_diff_corr_I, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr_I, bias_modis_diff_corr_I, bias_vpm_diff_corr_I, bias_bess_v1_diff_corr_I, bias_pmodel_diff_corr_I )

  tmp2 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr_II, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV )

  tmp3 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_III, bias_vpm_diff_corr_III, bias_bess_v1_diff_corr_III, bias_pmodel_diff_corr_III, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr_III, bias_modis_diff_corr_III, bias_vpm_diff_corr_III, bias_bess_v1_diff_corr_III, bias_pmodel_diff_corr_III )

  tmp4 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr_IV, bias_modis_diff_corr_IV, bias_vpm_diff_corr_IV, bias_bess_v1_diff_corr_IV, bias_pmodel_diff_corr_IV )

  tmp5 <- df_dday_8d_agg_norm %>% select( mysitename, date, fvar, infbin, insbin, bias_modis_diff_corr_V, bias_vpm_diff_corr_V, bias_bess_v1_diff_corr_V, bias_pmodel_diff_corr_V, gpp_obs ) %>% 
                                  gather( model, bias_diff_corr_V, bias_modis_diff_corr_V, bias_vpm_diff_corr_V, bias_bess_v1_diff_corr_V, bias_pmodel_diff_corr_V )

  out <- list( df_dday_8d_agg=df_dday_8d_agg, df_dday_8d_agg_norm=df_dday_8d_agg_norm, tmp=tmp, tmp0=tmp0, tmp1=tmp1, tmp2=tmp2, tmp3=tmp3, tmp4=tmp4, tmp5=tmp5 )                                  
  return( out )

}

get_numbers_biasreduction <- function( out ){
  ##------------------------------------------------
  ## Get numbers for paper: bias reduction
  ##------------------------------------------------
  bybin_orig <- out$tmp  %>% group_by( infbin ) %>% summarise( bias = mean(bias_diff, na.rm=TRUE), rmse = rmse( bias_diff + gpp_obs, gpp_obs) )
  bybin_corr <- out$tmp0 %>% group_by( infbin ) %>% summarise( bias = mean(bias_diff_corr, na.rm=TRUE), rmse = rmse( bias_diff_corr + gpp_obs, gpp_obs) )
  bybin_cor4 <- out$tmp4 %>% group_by( infbin ) %>% summarise( bias = mean(bias_diff_corr_IV, na.rm=TRUE), rmse = rmse( bias_diff_corr_IV + gpp_obs, gpp_obs) )

  print(paste("With fLUE, mean bias in lower four fLUE bins is reduced from", format(mean(bybin_orig$bias[1:4]), digits=3), "to", format(mean(bybin_corr$bias[1:4]), digits=3)))
  print(paste("With fLUE, RMSE in lower four fLUE bins is reduced from", format(mean(bybin_orig$rmse[1:4]), digits=3), "to", format(mean(bybin_corr$rmse[1:4]), digits=3)))

  print(paste("With IV, mean bias in lower four fLUE bins is reduced from", format(mean(bybin_orig$bias[1:4]), digits=3), "to", format(mean(bybin_cor4$bias[1:4]), digits=3)))
  print(paste("With IV, RMSE in lower four fLUE bins is reduced from", format(mean(bybin_orig$rmse[1:4]), digits=3), "to", format(mean(bybin_cor4$rmse[1:4]), digits=3)))
}

##------------------------------------------------
## Final plot, pooled models, corrected by fLUE
## showing:
## - uncorrected, normalised, pooled models
## - corrected by I, IV, and III, normalised, pooled models
## - corrected by fLUE, normalised, pooled models
##------------------------------------------------
xlim <- c(0.5,5.5)
ylim <- c(-6.5,6.5)

makepdf <- TRUE

plot_bias_resolved_flue <- function( tmp, tmp0, filn=NA ){
  ## pooled
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (mod.-obs., gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_diff ~ infbin, data = tmp, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_diff_corr ~ infbin, data = tmp0, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("pooled models, normalised", "pooled models, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_flue( tmp, tmp0, filn="fig/bias_resolved_fLUE.pdf")

plot_bias_resolved <- function( tmp, tmp0, tmp1, tmp4, tmp3, filn=NA, cex=1.0 ){
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1, mar=c(4,4,2,1) )
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (mod.-obs., gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_diff ~ infbin, data = tmp, at=1:5-0.2, col="tomato", boxwex=0.2, add=TRUE )
    myboxplot( bias_diff_corr ~ infbin, data = tmp0, at=1:5+0.0, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.2 )
    myboxplot( bias_diff_corr_I ~ infbin, data = tmp1, at=1:5+0.15, add=TRUE, col="springgreen1", axes=FALSE, boxwex=0.1 )
    myboxplot( bias_diff_corr_IV ~ infbin, data = tmp4, at=1:5+0.25, add=TRUE, col="springgreen3", axes=FALSE, boxwex=0.1 )
    myboxplot( bias_diff_corr_III ~ infbin, data = tmp3, at=1:5+0.35, add=TRUE, col="springgreen4", axes=FALSE, boxwex=0.1 )
    axis( 2, lwd=1.5 )
    axis( 4, lwd=1.5, labels=FALSE )
    box( lwd=1.5 )
    abline( h=0, lty=3 )
    legend("bottomleft", 
      c("pooled models, normalised", 
        "pooled models, normalised, corrected by fLUE", 
        "pooled models, normalised, corrected by I", 
        "pooled models, normalised, corrected by IV", 
        "pooled models, normalised, corrected by III"), 
      fill=c("tomato", "royalblue3", "springgreen1", "springgreen3", "springgreen4"), bty="n", cex=cex )
  if (!is.na(filn)) dev.off()  
}
# plot_bias_resolved( tmp, tmp0, tmp1, tmp4, tmp3, filn="fig/bias_resolved.pdf" )

plot_bias_resolved_pmodel <- function( df_dday_8d_agg_norm, filn=NA ){
  ## pmodel
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_pmodel_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_pmodel_diff_corr ~ infbin, data = df_dday_8d_agg_norm, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("P-model, normalised", "P-model, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_pmodel( df_dday_8d_agg_norm, filn="fig/bias_pmodel_resolved_fLUE.pdf" )

plot_bias_resolved_modis <- function( df_dday_8d_agg_norm, filn=NA ){
  ## modis
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_modis_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_modis_diff_corr ~ infbin, data = df_dday_8d_agg_norm, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("MODIS, normalised", "MODIS, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_modis( df_dday_8d_agg_norm, filn="fig/bias_modis_resolved_fLUE.pdf" )


plot_bias_resolved_vpm <- function( df_dday_8d_agg_norm, filn=NA ){
  ## vpm
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_vpm_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_vpm_diff_corr ~ infbin, data = df_dday_8d_agg_norm, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("VPM, normalised", "VPM, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_vpm( df_dday_8d_agg_norm, filn="fig/bias_vpm_resolved_fLUE.pdf" )

plot_bias_resolved_bess_v1 <- function( df_dday_8d_agg_norm, filn=NA ){
  ## bess
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_bess_v1_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_bess_v1_diff_corr ~ infbin, data = df_dday_8d_agg_norm, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("BESS, normalised", "BESS, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_bess_v1( df_dday_8d_agg_norm, filn="fig/bias_bess_v1_resolved_fLUE.pdf" )

plot_bias_resolved_mte <- function( df_dday_8d_agg_norm, filn=NA ){
  ## mte
  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)
    plot( xlim, ylim, type="n", ylim=ylim, xlim=xlim, xlab = "fLUE bin", ylab = expression( paste("bias (gC m"^-2, "d"^-1, ")" ) ), axes=FALSE )
    rect( 1:5-0.5, rep(ylim[1], 6), 1:5+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )
    myboxplot( bias_mte_diff ~ infbin, data = df_dday_8d_agg_norm, at=1:5-0.15, col="tomato", boxwex=0.3, add=TRUE )
    myboxplot( bias_mte_diff_corr ~ infbin, data = df_dday_8d_agg_norm, at=1:5+0.15, add=TRUE, col="royalblue3", axes=FALSE, boxwex=0.3 )
    abline( h=0, lty=3 )
    legend("bottomleft", c("MTE, normalised", "MTE, normalised, corrected by fLUE"), fill=c("tomato", "royalblue3"), bty="n")
  if (!is.na(filn)) dev.off()
}
# plot_bias_resolved_mte( df_dday_8d_agg_norm, filn="fig/bias_mte_resolved_fLUE.pdf" )



