syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

require(dplyr)
require(tidyr)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )
# source( "compl_df_flue_est.R" )

source("stress_quad_1sided.R")

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

##------------------------------------------------
## Manual
##------------------------------------------------
addboxes = FALSE 
##------------------------------------------------

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_8d_agg', 'df_dday_mte_agg', 'df_dday_bess_agg', 'df_dday_vpm_agg'

# ## Estimate soil moisture correction (adds column 'flue_est' to dataframe)
# load( "data/linearfit2_ratio.Rdata" )
# df_dday_agg       <- compl_df_flue_est( df_dday_agg   , linearfit2, x0_fix=0.9  )
# df_dday_8d_agg    <- compl_df_flue_est( df_dday_8d_agg, linearfit2, x0_fix=0.9  )

## define fLUE bins
nbins <- 11
max <- 1.1
binwidth <- max/nbins
fvarbins <- seq( from=0, to=max, by=binwidth )
xvals <- fvarbins[1:nbins]+binwidth/2

# ## define soil moisture bins
# nbins <- 10
# max <- 1.0
# binwidth <- max/nbins
# soilmbins <- seq( from=0, to=max, by=binwidth )
# xvals_soilm <- soilmbins[1:nbins]+binwidth/2

## bin data
df_dday_agg       <- df_dday_agg    %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
df_dday_8d_agg    <- df_dday_8d_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) ) %>%
                                        mutate( ifelse( is.nan(ratio_obs_mod_pmodel), NA, ratio_obs_mod_pmodel ) )

## load nice_agg to get data outside droughts
load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads nice_agg
load( "data/nice_nn_8d_agg_lue_obs_evi.Rdata" )   # loads nice_8d_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )

##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
## Merge mean annual alpha (AET/PET) values into this dataframe
load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
df_dday_agg    <- df_dday_agg    %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )
df_dday_8d_agg <- df_dday_8d_agg %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

df_dday_agg <- df_dday_agg %>% mutate(  ratio_obs_mod_pmodel_corr  = ifelse( flue_est_2>0, ratio_obs_mod_pmodel  / flue_est_2, NA ),
                                        ratio_obs_mod_bess_v1_corr = ifelse( flue_est_2>0, ratio_obs_mod_bess_v1 / flue_est_2, NA ),
                                        ratio_obs_mod_bess_v2_corr = ifelse( flue_est_2>0, ratio_obs_mod_bess_v2 / flue_est_2, NA )
                                      )
 
df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  ratio_obs_mod_modis_corr = ifelse( flue_est_2>0, ratio_obs_mod_modis / flue_est_2, NA ),
                                              ratio_obs_mod_vpm_corr   = ifelse( flue_est_2>0, ratio_obs_mod_vpm / flue_est_2, NA ) 
                                            )

## test 
for (sitename in unique(df_dday_agg$mysitename)){
  
  ## point cloud
  xlim <- c(0,1.1)
  ylim <- c(0,2)
  with( 
    filter( df_dday_agg, mysitename == sitename ),  # necessary to get useful bins with plot()
    plot( 
      fvar, 
      ratio_obs_mod_pmodel_corr, 
      xlab="",
      ylab="GPP observed / GPP modelled",
      xlim=xlim,
      ylim=ylim,
      main="",
      pch=16
    )
  )
  
  abline( h=1.0, lwd=0.5, lty=2 )
  abline( v=1.0, lwd=0.5, lty=2 )
  lines( c(-99,99), c(-99,99), col='black' )
  mtext( paste("P-model ", sitename), line=0.5, adj=0, font=2, cex=0.8 )
    
}


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

# pdf( "fig/bias_vs_fvar_boxes.pdf", width=sum(widths), height=sum(heights) )
# 
#   panel <- layout(
#                   order,
#                   widths=widths,
#                   heights=heights,
#                   TRUE
#                   )
#   # layout.show(panel)

  #---------------------------------------------------------
  # P-model
  #---------------------------------------------------------
    ## point cloud
    par( las=1, mar=c(2,4.5,2.5,rightmar) )
    xlim <- c(0,1.1)
    ylim <- c(0,2)
    with( 
          filter( df_dday_agg, ratio_obs_mod_pmodel_corr<5 ),  # necessary to get useful bins with plot()
          plot( 
                fvar, 
                ratio_obs_mod_pmodel_corr, 
                xlab="",
                ylab="GPP observed / GPP modelled",
                xlim=xlim,
                ylim=ylim,
                main="",
                type="n"
              )
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "P-model", line=0.5, adj=0, font=2, cex=0.8 )

    ## get distribution of bias within bins
    df_bins <- df_dday_agg %>%  group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_pmodel_corr) ) %>% 
                                summarise(  peak_corr      = getpeak(      ratio_obs_mod_pmodel_corr ), 
                                            uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_pmodel_corr, lev=0.75 ), 
                                            lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_pmodel_corr, lev=0.75 ),
                                            peak           = getpeak(      ratio_obs_mod_pmodel      ), 
                                            uhalfpeak      = getuhalfpeak( ratio_obs_mod_pmodel     , lev=0.75 ), 
                                            lhalfpeak      = getlhalfpeak( ratio_obs_mod_pmodel     , lev=0.75 )
                                             ) %>%
                                mutate( mids=xvals )

    ## plot uncorrected  
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

    ## plot CORRECTED distribution within bins
    rect( xvals-0.02, df_bins$lhalfpeak_corr, xvals+0.02, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
    with( df_bins, points( xvals, peak_corr, pch='-', col="springgreen3", cex=2 ) )
    
    # write stats into plot
    sub <- filter( df_dday_agg, is_drought_byvar==1 & ratio_obs_mod_pmodel_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_pmodel * gpp_pmodel )
    stats_s0 <- analyse_modobs( sub$gpp_pmodel, sub$gpp_obs, do.plot=FALSE )

    sub <- filter( df_dday_agg, is_drought_byvar==1 & ratio_obs_mod_pmodel_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_pmodel * gpp_pmodel )
    stats_s1 <- analyse_modobs( sub$gpp_pmodel * sub$flue_est_2, sub$gpp_obs, do.plot=FALSE )
    
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]
    mtext( paste0( "RMSE, s0 = ", format( stats_s0$rmse, digits = 3 ) ), adj=0.05, line=-1, cex=0.6 )
    mtext( paste0( "RMSE, s1 = ", format( stats_s1$rmse, digits = 3 ), " (-", format( (1 - stats_s1$rmse/stats_s0$rmse)*100 , digits = 3 ) ,"%)" ), adj=0.05, line=-2, cex=0.6 )
    
  #---------------------------------------------------------
  # MODIS
  #---------------------------------------------------------
    par( las=1, mar=c(2,2,2.5,rightmar) )
    xlim <- c(0,1.1)
    ylim <- c(0,2)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_modis_corr<5 ),  # necessary to get useful bins with plot()
          plot( 
                fvar, 
                ratio_obs_mod_modis_corr, 
                xlab="",
                ylab="",
                xlim=xlim,
                ylim=ylim,
                main="",
                type="n"
              ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "MOD17A2H", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    df_bins <- df_dday_8d_agg %>%  group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_modis_corr) ) %>% 
                                      summarise(  peak_corr      = getpeak(      ratio_obs_mod_modis_corr ), 
                                                  uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_modis_corr, lev=0.75 ), 
                                                  lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_modis_corr, lev=0.75 ),
                                                  peak           = getpeak(      ratio_obs_mod_modis      ), 
                                                  uhalfpeak      = getuhalfpeak( ratio_obs_mod_modis     , lev=0.75 ), 
                                                  lhalfpeak      = getlhalfpeak( ratio_obs_mod_modis     , lev=0.75 )
                                                   ) %>%
                                      mutate( mids=xvals )

    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

    ## plot CORRECTED distribution within bins
    rect( xvals-0.02, df_bins$lhalfpeak_corr, xvals+0.02, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
    with( df_bins, points( xvals, peak_corr, pch='-', col="springgreen3", cex=2 ) )

    # write stats into plot
    sub <- filter( df_dday_8d_agg, is_drought_byvar==1 & ratio_obs_mod_modis_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_modis * gpp_modis )
    stats_s0 <- analyse_modobs( sub$gpp_modis, sub$gpp_obs, do.plot=FALSE )

    sub <- filter( df_dday_8d_agg, is_drought_byvar==1 & ratio_obs_mod_modis_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_modis * gpp_modis )
    stats_s1 <- analyse_modobs( sub$gpp_modis * sub$flue_est_2, sub$gpp_obs, do.plot=FALSE )
    
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]
    mtext( paste0( "RMSE, s0 = ", format( stats_s0$rmse, digits = 3 ) ), adj=0.05, line=-1, cex=0.6 )
    mtext( paste0( "RMSE, s1 = ", format( stats_s1$rmse, digits = 3 ), " (-", format( (1 - stats_s1$rmse/stats_s0$rmse)*100 , digits = 3 ) ,"%)" ), adj=0.05, line=-2, cex=0.6 )
    
  #---------------------------------------------------------
  # BESS v1
  #---------------------------------------------------------
    par( las=1, mar=c(2,4.5,2.5,rightmar) )
    xlim <- c(0,1.1)
    ylim <- c(0,2)
    with( 
          filter( df_dday_agg, ratio_obs_mod_bess_v1_corr<5 ),  # necessary to get useful bins with plot()
          plot( 
                fvar, 
                ratio_obs_mod_bess_v1_corr, 
                xlab="",
                ylab="GPP observed / GPP modelled",
                xlim=xlim,
                ylim=ylim,
                main="",
                type="n"
              ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "BESS v1", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    df_dday_agg <- df_dday_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_agg %>% group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_bess_v1_corr) ) %>% 
                                    summarise(  peak_corr      = getpeak(      ratio_obs_mod_bess_v1_corr ), 
                                                uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_bess_v1_corr, lev=0.75 ), 
                                                lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_bess_v1_corr, lev=0.75 ),
                                                peak           = getpeak(      ratio_obs_mod_bess_v1      ), 
                                                uhalfpeak      = getuhalfpeak( ratio_obs_mod_bess_v1     , lev=0.75 ), 
                                                lhalfpeak      = getlhalfpeak( ratio_obs_mod_bess_v1     , lev=0.75 )
                                                 ) %>%
                                    mutate( mids=xvals )

    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
  
    ## plot CORRECTED distribution within bins
    rect( xvals-0.02, df_bins$lhalfpeak_corr, xvals+0.02, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
    with( df_bins, points( xvals, peak_corr, pch='-', col="springgreen3", cex=2 ) )

    # write stats into plot
    sub <- filter( df_dday_agg, is_drought_byvar==1 & ratio_obs_mod_bess_v1_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_bess_v1 * gpp_bess_v1 )
    stats_s0 <- analyse_modobs( sub$gpp_bess_v1, sub$gpp_obs, do.plot=FALSE )

    sub <- filter( df_dday_agg, is_drought_byvar==1 & ratio_obs_mod_bess_v1_corr<5 ) %>% mutate( gpp_obs = ratio_obs_mod_bess_v1 * gpp_bess_v1 )
    stats_s1 <- analyse_modobs( sub$gpp_bess_v1 * sub$flue_est_2, sub$gpp_obs, do.plot=FALSE )
    
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]
    mtext( paste0( "RMSE, s0 = ", format( stats_s0$rmse, digits = 3 ) ), adj=0.05, line=-1, cex=0.6 )
    mtext( paste0( "RMSE, s1 = ", format( stats_s1$rmse, digits = 3 ), " (-", format( (1 - stats_s1$rmse/stats_s0$rmse)*100 , digits = 3 ) ,"%)" ), adj=0.05, line=-2, cex=0.6 )
    
  #---------------------------------------------------------
  # BESS v2
  #---------------------------------------------------------
    par( las=1, mar=c(2,2,2.5,rightmar) )
    xlim <- c(0,1.1)
    ylim <- c(0,2)
    with( 
          filter( df_dday_agg, ratio_obs_mod_bess_v2<5 ),  # necessary to get useful bins with plot()
          plot( 
                fvar, 
                ratio_obs_mod_bess_v2, 
                xlab="",
                ylab="",
                xlim=xlim,
                ylim=ylim,
                main="",
                type="n"
              ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "BESS v2", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    df_dday_agg <- df_dday_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_agg %>%  group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_bess_v2) ) %>% 
                                summarise(  peak_corr      = getpeak(      ratio_obs_mod_bess_v2_corr ), 
                                            uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_bess_v2_corr, lev=0.75 ), 
                                            lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_bess_v2_corr, lev=0.75 ),
                                            peak           = getpeak(      ratio_obs_mod_bess_v2 ), 
                                            uhalfpeak      = getuhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 ), 
                                            lhalfpeak      = getlhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 )
                                             ) %>%
                                mutate( mids=xvals )
    
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
  
    ## plot CORRECTED distribution within bins
    rect( xvals-0.02, df_bins$lhalfpeak_corr, xvals+0.02, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
    with( df_bins, points( xvals, peak_corr, pch='-', col="springgreen3", cex=2 ) )

  #---------------------------------------------------------
  # VPM
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,rightmar) )
    xlim <- c(0,1.1)
    ylim <- c(0,2)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_vpm<5 ),  # necessary to get useful bins with plot()
          plot( 
                fvar, 
                ratio_obs_mod_vpm_corr, 
                xlab="fLUE",
                ylab="GPP observed / GPP modelled",
                xlim=xlim,
                ylim=ylim,
                main="",
                type="n"
              ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "VPM", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    df_dday_8d_agg <- df_dday_8d_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_8d_agg   %>% group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_vpm_corr) ) %>% 
                                    summarise(  peak_corr      = getpeak(      ratio_obs_mod_vpm_corr ), 
                                                uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_vpm_corr, lev=0.75 ), 
                                                lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_vpm_corr, lev=0.75 ),
                                                peak           = getpeak(      ratio_obs_mod_vpm ), 
                                                uhalfpeak      = getuhalfpeak( ratio_obs_mod_vpm, lev=0.75 ), 
                                                lhalfpeak      = getlhalfpeak( ratio_obs_mod_vpm, lev=0.75 )
                                                 ) %>%
                                    mutate( mids=xvals )

    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
  
    ## plot CORRECTED distribution within bins
    rect( xvals-0.02, df_bins$lhalfpeak_corr, xvals+0.02, df_bins$uhalfpeak_corr, col = add_alpha("springgreen3", 0.5), border = NA )
    with( df_bins, points( xvals, peak_corr, pch='-', col="springgreen3", cex=2 ) )

  #---------------------------------------------------------
  # MTE
  #---------------------------------------------------------
    par( las=1, mar=c(4,2,2.5,rightmar) )
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with plot()
          plot( 
                      fvar, 
                      ratio_obs_mod_mte, 
                      xlab="fLUE",
                      ylab="",
                      xlim=c(0,1.2),
                      ylim=ylim,
                      main="",
                      type="n"
                    ) 

        )
    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "FLUXCOM MTE", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    df_dday_8d_agg <- df_dday_8d_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_8d_agg %>% group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_mte) ) %>% 
                                  summarise(  peak_corr = getpeak(ratio_obs_mod_mte), 
                                              uhalfpeak_corr = getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), 
                                              lhalfpeak_corr = getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ),
                                              peak = getpeak(ratio_obs_mod_mte), 
                                              uhalfpeak = getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), 
                                              lhalfpeak = getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ),
                                              q25=quantile(ratio_obs_mod_mte, probs=0.25), 
                                              q75=quantile(ratio_obs_mod_mte, probs=0.75) 
                                              ) %>%
                                  complete( inbin, fill = list( peak=NA, uhalfpeak=NA, lhalfpeak=NA, q25=NA, q75=NA ) ) %>%
                                  mutate( mids=xvals )
    
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

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

  #   df_dday_8d_agg <- df_dday_8d_agg %>% mutate( inbin = cut( fvar, breaks = fvarbins ) )
    
  #   xvals <- fvarbins[1:nbins]+binwidth/2
  #   df_bins <- df_dday_8d_agg %>% group_by( inbin ) %>% filter( !is.na(inbin) & !is.na(ratio_obs_mod_rf) ) %>% 
  #     summarise( peak=getpeak(ratio_obs_mod_rf), uhalfpeak=getuhalfpeak( ratio_obs_mod_rf, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_rf, lev=0.75 ),
  #                q25=quantile(ratio_obs_mod_rf, probs=0.25), q75=quantile(ratio_obs_mod_rf, probs=0.75) ) %>%
  #     mutate( mids=xvals )
    
  #   rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("red", 0.5), border = NA )
  #   with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

  #   #---------------------------------------------------------
  #   # Distribution 
  #   #---------------------------------------------------------
  #   par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

  #   boxplot( filter( df_dday_8d_agg, dday < 0 )$ratio_obs_mod_rf, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
  #   abline( h=1.0, lwd=0.5, lty=2 )

# dev.off()




