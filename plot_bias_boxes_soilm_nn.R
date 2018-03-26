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

myboxplot <- function( ... ){
  bp <- boxplot( ..., staplewex = 0.0, whisklty = 1 )
  return(bp)
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

df_binned <- df_dday_agg %>% group_by( insoilmbin ) %>%
                             summarise( rmse_s0 = rmse( gpp_pmodel, gpp_obs ), rmse_s1 = rmse( gpp_pmodel * flue_est_2, gpp_obs ),
                                        meanbias_bess_v1 = median(bias_bess_v1, na.rm=TRUE),
                                        meanbias_bess_v2 = median(bias_bess_v2, na.rm=TRUE),
                                        meanbias_pmodel = median(bias_pmodel, na.rm=TRUE)
                                       ) %>%
                             mutate( reduction = (1 - rmse_s1 / rmse_s0 ) * 100 )

df_8d_binned <- df_dday_8d_agg %>% group_by( insoilmbin ) %>%
                                   summarise( meanbias_modis = median(bias_modis, na.rm=TRUE),
                                              meanbias_vpm = median(bias_vpm, na.rm=TRUE),
                                              meanbias_mte = median(bias_mte, na.rm=TRUE)
                                             )

##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
df_dday_agg <- df_dday_agg %>% mutate(  normdbias_bin0_pmodel = bias_pmodel / df_binned$meanbias_pmodel[nrow(df_binned)],
                                        normdbias_bin0_bess_v1 = bias_bess_v1 / df_binned$meanbias_bess_v1[nrow(df_binned)],
                                        normdbias_bin0_bess_v2 = bias_bess_v2 / df_binned$meanbias_bess_v2[nrow(df_binned)]
                                      )
 
df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  normdbias_bin0_modis = bias_modis / df_8d_binned$meanbias_modis[nrow(df_binned)],
                                              normdbias_bin0_vpm = bias_vpm / df_8d_binned$meanbias_vpm[nrow(df_binned)],
                                              normdbias_bin0_mte = bias_mte / df_8d_binned$meanbias_mte[nrow(df_binned)]
                                            )

##------------------------------------------------
## filter out some sites 
##------------------------------------------------
df_dday_agg    <- df_dday_agg    %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp", "AU-Fog", "AU-DaP", "AU-ASM", "IT-SRo", "US-SRG", "US-SRM") ) )
df_dday_8d_agg <- df_dday_8d_agg %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp") ) )

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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_pmodel ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_pmodel * flue_est_2 ) ~ insoilmbin, 
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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_modis ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_modis * flue_est_2 ) ~ insoilmbin, 
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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_bess_v1 ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_bess_v1 * flue_est_2 ) ~ insoilmbin, 
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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_bess_v2 ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_bess_v2 * flue_est_2 ) ~ insoilmbin, 
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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_vpm ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_vpm * flue_est_2 ) ~ insoilmbin, 
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
    bp1 <- myboxplot( 
                    ( normdbias_bin0_mte ) ~ insoilmbin, 
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
                    ylab="GPP modelled / observed, norm.", 
                    xlab="soil moisture bins"  
                    )
    bp4 <- myboxplot( 
                    ( normdbias_bin0_mte * flue_est_2 ) ~ insoilmbin, 
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

  #   myboxplot( filter( df_dday_8d_agg, dday < 0 )$ratio_obs_mod_rf, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
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
  bp1 <- myboxplot( 
                  ( normdbias_bin0_pmodel ) ~ insoilmbin, 
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
                  ylab="GPP modelled / observed, norm.", 
                  xlab="soil moisture bins",
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_1 ) ~ insoilmbin, 
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
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_2 ) ~ insoilmbin, 
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
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_3 ) ~ insoilmbin, 
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
  bp1 <- myboxplot( 
                  ( normdbias_bin0_pmodel ) ~ infvarbin, 
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
                  ylab="GPP modelled / observed, norm.", 
                  xlab="fLUE bins",
                  staplewex = 0.0,
                  whisklty = 1
                  )
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_1 ) ~ infvarbin, 
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
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_2 ) ~ infvarbin, 
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
  bp4 <- myboxplot( 
                  ( normdbias_bin0_pmodel * flue_est_3 ) ~ infvarbin, 
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

