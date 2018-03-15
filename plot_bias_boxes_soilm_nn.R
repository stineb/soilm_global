syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

require(dplyr)
require(LSD)
require(tidyr)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )
source( "compl_df_flue_est.R" )

getpeak <- function( vec ) {
  vec  <- vec[!is.na(vec)]
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[ dens$y==max(dens$y) ]
  return( peak )   
}

getlhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  vec  <- vec[!is.na(vec)]
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x<peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
} 

getuhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  vec  <- vec[!is.na(vec)]
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x>peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
}

rmse <- function( mod, obs ){
  rmse <- sqrt( mean( (mod - obs)^2, na.rm = TRUE ) )
  return(rmse)
}

##------------------------------------------------
## Manual
##------------------------------------------------
addboxes = FALSE 
##------------------------------------------------

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_8d_agg', 'df_dday_mte_agg', 'df_dday_bess_agg', 'df_dday_vpm_agg'

## Use only sites where NN method worked (i.e. that had clear and identifiable soil moisture limitation)
successcodes <- read.csv( paste0( myhome, "/sofun/utils_sofun/analysis_sofun/fluxnet2015/successcodes.csv" ), as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename
df_dday_agg <- df_dday_agg %>% filter( mysitename %in% do.sites )

## define soil moisture bins
nbins <- 4
max <- 1.0
binwidth <- max/nbins
soilmbins <- seq( from=0, to=max, by=binwidth )
xvals <- soilmbins[1:nbins]+binwidth/2

## bin data
df_dday_agg       <- df_dday_agg    %>% mutate( insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ), infvarbin = cut( as.numeric(fvar), breaks = soilmbins ), gpp_obs = gpp_pmodel / bias_pmodel ) 
df_dday_8d_agg    <- df_dday_8d_agg %>% mutate( insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) %>%
                                        mutate( ifelse( is.nan(ratio_obs_mod_pmodel), NA, ratio_obs_mod_pmodel ) )

df_rmse <- df_dday_agg %>% group_by( insoilmbin ) %>%
                           summarise( rmse_s0 = rmse( gpp_pmodel, gpp_obs ), rmse_s1 = rmse( gpp_pmodel * flue_est_2, gpp_obs ) ) %>%
                           mutate( reduction = (1 - rmse_s1 / rmse_s0 ) * 100 )

##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
df_dday_agg <- df_dday_agg %>% mutate(  ratio_obs_mod_pmodel_corr = ratio_obs_mod_pmodel / flue_est_2,
                                        ratio_obs_mod_bess_v1_corr = ratio_obs_mod_bess_v1 / flue_est_2,
                                        ratio_obs_mod_bess_v2_corr = ratio_obs_mod_bess_v2 / flue_est_2
                                      )
 
df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  ratio_obs_mod_modis_corr = ratio_obs_mod_modis / flue_est_2,
                                              ratio_obs_mod_vpm_corr = ratio_obs_mod_vpm / flue_est_2 
                                            )

## quick and dirty
# bp1 <- boxplot( log( bias_pmodel ) ~ insoilmbin, data=ddf, col="tomato", las=1, outline = FALSE, na.rm=TRUE, add=FALSE, at=(soilmbins[1:nbins]+1/nbins*0.35), boxwex=0.052, xlim=c(0,1), ylab="log of obs./mod.", xlab="soil moisture bins"  )
# bp4 <- boxplot( log( bias_pmodel * flue_est ) ~ insoilmbin, data=ddf, col="springgreen3", las=1, outline = FALSE, na.rm=TRUE, add=TRUE, at=(soilmbins[1:nbins]+1/nbins*0.65), boxwex=0.052, axes=FALSE, xlim=c(0,1) )
# abline( h=0.0, lty=3 )
# legend("bottomleft", c("P-model", "P-model, corrected"), bty="n", fill=c("tomato", "springgreen3") )

# ## get distribution of bias within bins
# df_bins <- df_dday_agg %>%  group_by( insoilmbin ) %>% filter( !is.na(insoilmbin) & !is.na(ratio_obs_mod_pmodel_corr) ) %>% 
#                             summarise(  peak_corr      = getpeak(      ratio_obs_mod_pmodel_corr ), 
#                                         uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_pmodel_corr, lev=0.75 ), 
#                                         lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_pmodel_corr, lev=0.75 ),
#                                         peak           = getpeak(      ratio_obs_mod_pmodel      ), 
#                                         uhalfpeak      = getuhalfpeak( ratio_obs_mod_pmodel     , lev=0.75 ), 
#                                         lhalfpeak      = getlhalfpeak( ratio_obs_mod_pmodel     , lev=0.75 )
#                                          ) %>%
#                             mutate( mids=xvals )

# ## plot uncorrected  
# rect( xvals-0.04, df_bins$lhalfpeak, xvals+0.00, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
# with( df_bins, points( xvals-0.02, peak, pch='-', col="red", cex=2 ) )

# ## plot CORRECTED distribution within bins
# rect( xvals-0.0, df_bins$lhalfpeak_corr, xvals+0.04, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
# with( df_bins, points( xvals+0.02, peak_corr, pch='-', col="springgreen3", cex=2 ) )

##------------------------------------------------
## filter out some sites 
##------------------------------------------------
# df_dday_agg    <- df_dday_agg    %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp") ) )
# df_dday_8d_agg <- df_dday_8d_agg %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp") ) )

##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------

## panel setup
magn <- 3
nrows <- 3
ncols <- 2
widths <- 0.9*c(magn, 0.9*magn )
rightmar <- 1
heights <- 1.2*c(0.6*magn,0.6*magn,0.68*magn)
order <- matrix(seq(ncols*nrows),nrows,ncols,byrow=TRUE)

pdf( "fig/bias_vs_fvar_boxes_soilm.pdf", width=sum(widths), height=sum(heights) )

  panel <- layout(
                  order,
                  widths=widths,
                  heights=heights,
                  TRUE
                  )
  # layout.show(panel)

  #---------------------------------------------------------
  # P-model
  #---------------------------------------------------------
    ## point cloud
    par( las=1, mar=c(2,4.5,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_pmodel ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.052, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="observed/modelled", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_pmodel * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.052, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "P-model", line=0.5, adj=0, font=2, cex=0.8 )
    legend("topright", c("s0", "s1b"), bty="n", fill=c("tomato", "springgreen3") )

  #---------------------------------------------------------
  # MODIS
  #---------------------------------------------------------
    par( las=1, mar=c(2,2,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_modis ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.052, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="log of obs./mod.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_modis * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.052, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "MOD17A2H", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # BESS v1
  #---------------------------------------------------------
    par( las=1, mar=c(2,4.5,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_bess_v1 ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.052, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="log of obs./mod.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_bess_v1 * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.052, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "BESS v1", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # BESS v2
  #---------------------------------------------------------
    par( las=1, mar=c(2,2,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_bess_v2 ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.052, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="log of obs./mod.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_bess_v2 * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.052, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "BESS v2", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # VPM
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_vpm ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.055, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="log of obs./mod.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_vpm * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.055, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "VPM", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # MTE
  #---------------------------------------------------------
    par( las=1, mar=c(4,2,2.5,rightmar) )
    xlim <- c(0,1)
    ylim <- c(0,6)
    bp1 <- boxplot( 
                    ( bias_mte ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="tomato", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=FALSE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.35), 
                    boxwex=0.055, 
                    xlim=xlim, 
                    ylim=ylim,
                    ylab="log of obs./mod.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- boxplot( 
                    ( bias_mte * flue_est_2 ) ~ insoilmbin, 
                    data=df_dday_8d_agg, 
                    col="springgreen3", 
                    las=1, 
                    outline = FALSE, 
                    na.rm=TRUE, 
                    add=TRUE, 
                    at=(soilmbins[1:nbins]+1/nbins*0.65), 
                    boxwex=0.055, 
                    axes=FALSE, 
                    xlim=xlim,
                    ylim=ylim 
                    )
    abline( h=1.0, lty=3 )
    mtext( "FLUXCOM MTE", line=0.5, adj=0, font=2, cex=0.8 )

  # #---------------------------------------------------------
  # # MTE-RF
  # #---------------------------------------------------------
  #   par( las=1, mar=c(4,4.5,2.5,0) )
  #   with( 
  #         filter( df_dday_8d_agg, ratio_obs_mod_rf<5 ),  # necessary to get useful bins with plot()
  #         plot( 
  #               fvar, 
  #               ratio_obs_mod_rf, 
  #               xlab="fLUE",
  #               ylab="GPP observed / GPP modelled",
  #               xlim=c(0,1.2),
  #               ylim=ylim,
  #               cex=1.2,
  #               pch=16,
  #               col=add_alpha("black", 0.3),
  #               main=""
  #             ) 

  #       )
  #   abline( h=1.0, lwd=0.5, lty=2 )
  #   abline( v=1.0, lwd=0.5, lty=2 )
  #   lines( c(-99,99), c(-99,99), col='black' )
  #   mtext( "FLUXCOM MTE-RF", line=1, adj=0.5 )

  #   df_dday_8d_agg <- df_dday_8d_agg %>% mutate( insoilmbin = cut( fvar, breaks = fvarbins ) )
    
  #   xvals <- fvarbins[1:nbins]+binwidth/2
  #   df_bins <- df_dday_8d_agg %>% group_by( insoilmbin ) %>% filter( !is.na(insoilmbin) & !is.na(ratio_obs_mod_rf) ) %>% 
  #     summarise( peak=getpeak(ratio_obs_mod_rf), uhalfpeak=getuhalfpeak( ratio_obs_mod_rf, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_rf, lev=0.75 ),
  #                q25=quantile(ratio_obs_mod_rf, probs=0.25), q75=quantile(ratio_obs_mod_rf, probs=0.75) ) %>%
  #     mutate( mids=xvals )
    
  #   rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
  #   with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

  #   #---------------------------------------------------------
  #   # Distribution 
  #   #---------------------------------------------------------
  #   par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

  #   boxplot( filter( df_dday_8d_agg, dday < 0 )$ratio_obs_mod_rf, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
  #   abline( h=1.0, lwd=0.5, lty=2 )

dev.off()


pdf("fig/bias_vs_fvar_boxes_soilm_pmodel.pdf", width = 6, height = 5 )
#---------------------------------------------------------
# P-model
#---------------------------------------------------------
  ## point cloud
  par( las=1, mar=c(4,4.5,2.5,rightmar) )
  xlim <- c(0,1)
  ylim <- c(0,6.3)
  bp1 <- boxplot( 
                  ( bias_pmodel ) ~ insoilmbin, 
                  data=df_dday_agg, 
                  col="tomato", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=FALSE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.35), 
                  boxwex=0.052, 
                  xlim=xlim, 
                  ylim=ylim,
                  ylab="GPP modelled / observed", 
                  xlab="soil moisture bins",
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_1 ) ~ insoilmbin, 
                  data=df_dday_agg, 
                  col="springgreen1", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.55), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_2 ) ~ insoilmbin, 
                  data=df_dday_agg, 
                  col="springgreen3", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.65), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_3 ) ~ insoilmbin, 
                  data=df_dday_agg, 
                  col="springgreen4", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.75), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  abline( h=1.0, lty=3 )
  legend("topright", c("s0", "s1a", "s1b", "s1c"), bty="n", fill=c("tomato", "springgreen1", "springgreen3", "springgreen4") )

dev.off()


pdf("fig/bias_vs_fvar_boxes_flue_pmodel.pdf", width = 6, height = 5 )
#---------------------------------------------------------
# P-model
#---------------------------------------------------------
  ## point cloud
  par( las=1, mar=c(4,4.5,2.5,rightmar) )
  xlim <- c(0,1)
  ylim <- c(0,20)
  bp1 <- boxplot( 
                  ( bias_pmodel ) ~ infvarbin, 
                  data=df_dday_agg, 
                  col="tomato", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=FALSE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.35), 
                  boxwex=0.052, 
                  xlim=xlim, 
                  ylim=ylim,
                  ylab="GPP modelled / observed", 
                  xlab="fLUE bins",
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_1 ) ~ infvarbin, 
                  data=df_dday_agg, 
                  col="springgreen1", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.55), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_2 ) ~ infvarbin, 
                  data=df_dday_agg, 
                  col="springgreen3", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.65), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- boxplot( 
                  ( bias_pmodel * flue_est_3 ) ~ infvarbin, 
                  data=df_dday_agg, 
                  col="springgreen4", 
                  las=1, 
                  outline = FALSE, 
                  na.rm=TRUE, 
                  add=TRUE, 
                  at=(soilmbins[1:nbins]+1/nbins*0.75), 
                  boxwex=0.02, 
                  axes=FALSE, 
                  xlim=xlim,
                  ylim=ylim,
                  staplewex = 0.0,
                  whisklty = 1
                  )
  abline( h=1.0, lty=3 )
  legend("topright", c("s0", "s1a", "s1b", "s1c"), bty="n", fill=c("tomato", "springgreen1", "springgreen3", "springgreen4") )

dev.off()
