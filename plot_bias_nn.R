syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(LSD, quietly = TRUE, warn.conflicts = FALSE)

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

## group data by fLUE
nbins <- 10
binwidth <- 1.0/nbins
fvarbins <- seq( from=0, to=1, by=binwidth )
xvals <- fvarbins[1:nbins]+binwidth/2

df_dday_agg       <- df_dday_agg    %>% mutate( infvarbin = cut( fvar_smooth, breaks = fvarbins ) )
df_dday_8d_agg    <- df_dday_8d_agg %>% mutate( infvarbin = cut( fvar_smooth, breaks = fvarbins ) ) %>%
                                        mutate( ifelse( is.nan(ratio_obs_mod_pmodel), NA, ratio_obs_mod_pmodel ) )

## load nice_agg to get data outside droughts
load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads nice_agg
load( "data/nice_nn_8d_agg_lue_obs_evi.Rdata" )   # loads nice_8d_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )

##------------------------------------------------
## correct ratio with estimated fLUE
##------------------------------------------------
# df_dday_agg <- df_dday_agg %>% mutate(  ratio_obs_mod_pmodel_corr = ratio_obs_mod_pmodel,
#                                         ratio_obs_mod_bess_v1_corr = ratio_obs_mod_bess_v1
#                                       )

# df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  ratio_obs_mod_modis_corr = ratio_obs_mod_modis,
#                                               ratio_obs_mod_vpm_corr = ratio_obs_mod_vpm 
#                                             )

## Merge mean annual alpha (AET/PET) values into this dataframe
# load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
# df_dday_agg    <- df_dday_agg    %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )
# df_dday_8d_agg <- df_dday_8d_agg %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

# df_dday_agg <- df_dday_agg %>% mutate(  flue_est_3 = stress_quad_1sided_alpha( soilm_mean, meanalpha, x0 = 0.9, apar = -0.5055405, bpar = 0.8109020 ) ) %>%
#                                mutate(  ratio_obs_mod_pmodel_corr = ratio_obs_mod_pmodel, # / flue_est_3,
#                                         ratio_obs_mod_bess_v1_corr = ratio_obs_mod_bess_v1, # / flue_est_3,
#                                         ratio_obs_mod_bess_v2_corr = ratio_obs_mod_bess_v2 # / flue_est_3
#                                       )
 
# df_dday_8d_agg <- df_dday_8d_agg %>% mutate(  flue_est_3 = stress_quad_1sided_alpha( soilm_mean, meanalpha, x0 = 0.9, apar = -0.5055405, bpar = 0.8109020 ) ) %>%
#                                      mutate(  ratio_obs_mod_modis_corr = ratio_obs_mod_modis, # / flue_est_3,
#                                               ratio_obs_mod_vpm_corr = ratio_obs_mod_vpm, # / flue_est_3,
#                                               ratio_obs_mod_mte_corr = ratio_obs_mod_mte # / flue_est_3 
#                                             )

##------------------------------------------------
## filter out some sites 
##------------------------------------------------
# df_dday_agg    <- df_dday_agg    %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp") ) )
# df_dday_8d_agg <- df_dday_8d_agg %>% filter( !( mysitename %in% c("US-Var", "IT-Noe", "FR-Pue", "AU-Stp") ) )

##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------

ylim <- c(0,3.5)

## panel setup
magn <- 3
nrows <- 2
if (addboxes){
  ncols <- 4
  widths <- 0.9*c(magn, 0.2*magn, 0.85*magn, 0.2*magn )
  rightmar <- 0
} else {
  ncols <- 2
  widths <- 0.9*c(magn, 0.85*magn )
  rightmar <- 1
}
heights <- 1.2*c(0.6*magn,0.68*magn)
order <- matrix(seq(ncols*nrows),nrows,ncols,byrow=TRUE)

pdf( "fig/bias_vs_fvar.pdf", width=sum(widths), height=sum(heights) )

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
    xlim <- c(0,1.2)
    with( 
          filter( df_dday_agg, ratio_obs_mod_pmodel<10 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar_smooth, 
                      ratio_obs_mod_pmodel, 
                      xlab="",
                      ylab="GPP observed/modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main=""
                    )
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "P-model", line=0.5, adj=0, font=2, cex=0.8 )

    ## get distribution of bias within bins
    if (addboxes){

      df_bins <- df_dday_agg %>%  group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_pmodel) ) %>% 
                                  summarise(  peak      = getpeak(      ratio_obs_mod_pmodel ), 
                                              uhalfpeak = getuhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ), 
                                              lhalfpeak = getlhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ) ) %>%
                                  mutate( mids=xvals )

      ## plot uncorrected  
      rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
      with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

      # ## plot CORRECTED distribution within bins
      # rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
      # with( df_bins, points( xvals, peak, pch='-', col="springgreen", cex=2 ) )

    }

    # ## use only data during droughts for stats
    # sub <- filter( df_dday_agg, is_drought_byvar==1 ) 
    # stats <- analyse_modobs( sub$ratio_obs_mod_pmodel, sub$fvar, do.plot=FALSE )

    # # write stats into plot
    # x0 <- 0.05*xlim[2]
    # y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    # text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    # text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )
    
    ## draw the legend
    legend( "topleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", cex=0.8 )

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    if (addboxes){
      vec <- df_dday_agg %>% filter( dday < 0 & !is.na(ratio_obs_mod_pmodel) ) %>% select( ratio_obs_mod_pmodel )
      dens <- density( vec$ratio_obs_mod_pmodel )
      df_before <- vec %>% 
                      summarise( 
                        peak=getpeak(ratio_obs_mod_pmodel), 
                        uhalfpeak=getuhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ), 
                        lhalfpeak=getlhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ),
                        q25=quantile(ratio_obs_mod_pmodel, probs=0.25), 
                        q75=quantile(ratio_obs_mod_pmodel, probs=0.75) 
                        )
      par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
      xlim <- c(0,1)
      plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
      rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
      with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
      lines( dens$y/max(dens$y), dens$x )
      abline( h=1.0, lty=3 )
    } 
    # else {
    #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
    #   xlim <- c(0,1)
    #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
    # }
    
  #---------------------------------------------------------
  # MODIS
  #---------------------------------------------------------
    par( las=1, mar=c(2,2,2.5,rightmar) )
    xlim <- c(0,1.2)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_modis<10 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar_smooth, 
                      ratio_obs_mod_modis, 
                      xlab="",
                      ylab="",
                      xlim=xlim,
                      ylim=ylim,
                      cexplot=1.2,
                      main=""
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "MOD17A2H", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    if (addboxes){
      df_bins <- df_dday_8d_agg %>%  group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_modis) ) %>% 
                                        summarise(  peak           = getpeak(      ratio_obs_mod_modis ), 
                                                    uhalfpeak      = getuhalfpeak( ratio_obs_mod_modis, lev=0.75 ), 
                                                    lhalfpeak      = getlhalfpeak( ratio_obs_mod_modis, lev=0.75 ) ) %>%
                                        mutate( mids=xvals )

      rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
      with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
    }

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    if (addboxes){
      vec <- df_dday_8d_agg %>% filter( dday < 0 & !is.nan(ratio_obs_mod_modis) & !is.na(ratio_obs_mod_modis) ) %>% select( ratio_obs_mod_modis )
      dens <- density(vec$ratio_obs_mod_modis)
      df_before <- vec %>% 
                      summarise( 
                        peak=getpeak(ratio_obs_mod_modis), 
                        uhalfpeak=getuhalfpeak( ratio_obs_mod_modis, lev=0.75 ), 
                        lhalfpeak=getlhalfpeak( ratio_obs_mod_modis, lev=0.75 ),
                        q25=quantile(ratio_obs_mod_modis, probs=0.25), 
                        q75=quantile(ratio_obs_mod_modis, probs=0.75) 
                        )
      par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
      xlim <- c(0,1)
      plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
      rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
      with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
      lines( dens$y/max(dens$y), dens$x )
      abline( h=1.0, lty=3 )
    } 
    # else {
    #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
    #   xlim <- c(0,1)
    #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
    # }


  #---------------------------------------------------------
  # BESS v1
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,rightmar) )
    xlim <- c(0,1.2)
    with( 
          filter( df_dday_agg, ratio_obs_mod_bess_v1<10 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar_smooth, 
                      ratio_obs_mod_bess_v1, 
                      xlab="",
                      ylab="GPP observed/modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main=""
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "BESS v1", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    if (addboxes){
      df_dday_agg <- df_dday_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
      xvals <- fvarbins[1:nbins]+binwidth/2
      df_bins <- df_dday_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_bess_v1) ) %>% 
                                      summarise(  peak           = getpeak(      ratio_obs_mod_bess_v1 ), 
                                                  uhalfpeak      = getuhalfpeak( ratio_obs_mod_bess_v1, lev=0.75 ), 
                                                  lhalfpeak      = getlhalfpeak( ratio_obs_mod_bess_v1, lev=0.75 ) ) %>%
                                      mutate( mids=xvals )

      rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
      with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
    }
    
    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    if (addboxes){
      vec <- df_dday_agg %>% filter( dday < 0 & !is.nan(ratio_obs_mod_bess_v1) & !is.na(ratio_obs_mod_bess_v1) ) %>% select( ratio_obs_mod_bess_v1 )
      dens <- density(vec$ratio_obs_mod_bess_v1)
      df_before <- vec %>% 
                      summarise( 
                        peak=getpeak(ratio_obs_mod_bess_v1), 
                        uhalfpeak=getuhalfpeak( ratio_obs_mod_bess_v1, lev=0.75 ), 
                        lhalfpeak=getlhalfpeak( ratio_obs_mod_bess_v1, lev=0.75 ),
                        q25=quantile(ratio_obs_mod_bess_v1, probs=0.25), 
                        q75=quantile(ratio_obs_mod_bess_v1, probs=0.75) 
                        )
      par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
      xlim <- c(0,1)
      plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
      rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
      with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
      lines( dens$y/max(dens$y), dens$x )
      abline( h=1.0, lty=3 )
    } 
    # else {
    #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
    #   xlim <- c(0,1)
    #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
    # }


  # #---------------------------------------------------------
  # # BESS v2
  # #---------------------------------------------------------
  #   par( las=1, mar=c(2,2,2.5,rightmar) )
  #   xlim <- c(0,1.2)
  #   with( 
  #         filter( df_dday_agg, ratio_obs_mod_bess_v2<10 ),  # necessary to get useful bins with heatscatter()
  #         heatscatter( 
  #                     fvar_smooth, 
  #                     ratio_obs_mod_bess_v2, 
  #                     xlab="",
  #                     ylab="",
  #                     xlim=xlim,
  #                     ylim=ylim,
  #                     main=""
  #                   ) 
  #       )

  #   abline( h=1.0, lwd=0.5, lty=2 )
  #   abline( v=1.0, lwd=0.5, lty=2 )
  #   lines( c(-99,99), c(-99,99), col='red' )
  #   mtext( "BESS v2", line=0.5, adj=0, font=2, cex=0.8 )

  #   ## add boxes for distribution within bins
  #   if (addboxes){
  #     df_dday_agg <- df_dday_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
  #     xvals <- fvarbins[1:nbins]+binwidth/2
  #     df_bins <- df_dday_agg %>%  group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_bess_v2) ) %>% 
  #                                 summarise(  peak           = getpeak(      ratio_obs_mod_bess_v2 ), 
  #                                             uhalfpeak      = getuhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 ), 
  #                                             lhalfpeak      = getlhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 ) ) %>%
  #                                 mutate( mids=xvals )
      
  #     rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
  #     with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
  #   }
    
  #   #---------------------------------------------------------
  #   # Distribution 
  #   #---------------------------------------------------------
  #   if (addboxes){
  #     vec <- df_dday_agg %>% filter( dday < 0 & !is.nan(ratio_obs_mod_bess_v2) & !is.na(ratio_obs_mod_bess_v2) ) %>% select( ratio_obs_mod_bess_v2 )
  #     dens <- density(vec$ratio_obs_mod_bess_v2)
  #     df_before <- vec %>% 
  #                     summarise( 
  #                       peak=getpeak(ratio_obs_mod_bess_v2), 
  #                       uhalfpeak=getuhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 ), 
  #                       lhalfpeak=getlhalfpeak( ratio_obs_mod_bess_v2, lev=0.75 ),
  #                       q25=quantile(ratio_obs_mod_bess_v2, probs=0.25), 
  #                       q75=quantile(ratio_obs_mod_bess_v2, probs=0.75) 
  #                       )
  #     par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
  #     xlim <- c(0,1)
  #     plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
  #     rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
  #     with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
  #     lines( dens$y/max(dens$y), dens$x )
  #     abline( h=1.0, lty=3 )
  #   } 
  #   # else {
  #   #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
  #   #   xlim <- c(0,1)
  #   #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
  #   # }


  #---------------------------------------------------------
  # VPM
  #---------------------------------------------------------
    par( las=1, mar=c(4,2,2.5,rightmar) )
    xlim <- c(0,1.2)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_vpm<10 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar_smooth, 
                      ratio_obs_mod_vpm, 
                      xlab="fLUE",
                      ylab="",
                      xlim=xlim,
                      ylim=ylim,
                      cexplot=1.2,
                      main=""
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "VPM", line=0.5, adj=0, font=2, cex=0.8 )

    ## add boxes for distribution within bins
    if (addboxes){
      df_dday_8d_agg <- df_dday_8d_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
      xvals <- fvarbins[1:nbins]+binwidth/2
      df_bins <- df_dday_8d_agg   %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_vpm) ) %>% 
                                      summarise(  peak           = getpeak(      ratio_obs_mod_vpm ), 
                                                  uhalfpeak      = getuhalfpeak( ratio_obs_mod_vpm, lev=0.75 ), 
                                                  lhalfpeak      = getlhalfpeak( ratio_obs_mod_vpm, lev=0.75 ) ) %>%
                                      mutate( mids=xvals )

      rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
      with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
    }
    
    #---------------------------------------------------------
    # Distribution during days before drought onset
    #---------------------------------------------------------
    if (addboxes){
      vec <- df_dday_8d_agg %>% filter( dday < 0 & !is.nan(ratio_obs_mod_vpm) & !is.na(ratio_obs_mod_vpm) ) %>% select( ratio_obs_mod_vpm )
      dens <- density(vec$ratio_obs_mod_vpm)
      df_before <- vec %>% 
                      summarise( 
                        peak=getpeak(ratio_obs_mod_vpm), 
                        uhalfpeak=getuhalfpeak( ratio_obs_mod_vpm, lev=0.75 ), 
                        lhalfpeak=getlhalfpeak( ratio_obs_mod_vpm, lev=0.75 ),
                        q25=quantile(ratio_obs_mod_vpm, probs=0.25), 
                        q75=quantile(ratio_obs_mod_vpm, probs=0.75) 
                        )
      par( las=1, mar=c(4,0.2,2.5,2), xpd=FALSE )
      xlim <- c(0,1)
      plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
      rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
      with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
      lines( dens$y/max(dens$y), dens$x )
      abline( h=1.0, lty=3 )
    } 
    # else {
    #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
    #   xlim <- c(0,1)
    #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
    # }

  # #---------------------------------------------------------
  # # MTE
  # #---------------------------------------------------------
  #   par( las=1, mar=c(4,2,2.5,rightmar) )
  #   with( 
  #         filter( df_dday_8d_agg, ratio_obs_mod_mte<10 ),  # necessary to get useful bins with heatscatter()
  #         heatscatter( 
  #                     fvar_smooth, 
  #                     ratio_obs_mod_mte, 
  #                     xlab="fLUE",
  #                     ylab="",
  #                     xlim=c(0,1.2),
  #                     ylim=ylim,
  #                     cexplot=1.2,

  #                     main=""
  #                   ) 

  #       )
  #   abline( h=1.0, lwd=0.5, lty=2 )
  #   abline( v=1.0, lwd=0.5, lty=2 )
  #   lines( c(-99,99), c(-99,99), col='red' )
  #   mtext( "FLUXCOM MTE", line=0.5, adj=0, font=2, cex=0.8 )

  #   ## add boxes for distribution within bins
  #   if (addboxes){
  #     df_dday_8d_agg <- df_dday_8d_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
  #     xvals <- fvarbins[1:nbins]+binwidth/2
  #     df_bins <- df_dday_8d_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_mte) ) %>% 
  #       summarise(  peak=getpeak(ratio_obs_mod_mte), uhalfpeak=getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ),
  #                   q25=quantile(ratio_obs_mod_mte, probs=0.25), q75=quantile(ratio_obs_mod_mte, probs=0.75) ) %>%
  #       complete( infvarbin, fill = list( peak=NA, uhalfpeak=NA, lhalfpeak=NA, q25=NA, q75=NA ) ) %>%
  #       mutate( mids=xvals )
      
  #     rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
  #     with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
  #   }    

  #   #---------------------------------------------------------
  #   # Distribution 
  #   #---------------------------------------------------------
  #   if (addboxes){
  #     vec <- df_dday_8d_agg %>% filter( dday < 0 & !is.na(ratio_obs_mod_mte) & !is.nan(ratio_obs_mod_mte) ) %>% select( ratio_obs_mod_mte )
  #     dens <- density(vec$ratio_obs_mod_mte)
  #     df_before <- vec %>% 
  #                     summarise( 
  #                       peak=getpeak(ratio_obs_mod_mte), 
  #                       uhalfpeak=getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), 
  #                       lhalfpeak=getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ),
  #                       q25=quantile(ratio_obs_mod_mte, probs=0.25), 
  #                       q75=quantile(ratio_obs_mod_mte, probs=0.75) 
  #                       )
  #     par( las=1, mar=c(4,0.2,2.5,2), xpd=FALSE )
  #     xlim <- c(0,1)
  #     plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )
  #     rect( 0, df_before$lhalfpeak, 0.75, df_before$uhalfpeak, col = 'grey50' )
  #     with( df_before, points( 0.375, peak, pch='-', col="red", cex=4 ) )
  #     lines( dens$y/max(dens$y), dens$x )
  #     abline( h=1.0, lty=3 )
  #   } 
  #   # else {
  #   #   par( las=1, mar=c(2,0.2,2.5,2), xpd=FALSE )
  #   #   xlim <- c(0,1)
  #   #   plot( xlim, ylim, type="n", axes=FALSE, xlab="", ylab="" )      
  #   # }

  # #---------------------------------------------------------
  # # MTE-RF
  # #---------------------------------------------------------
  #   par( las=1, mar=c(4,4.5,2.5,0) )
  #   with( 
  #         filter( df_dday_8d_agg, ratio_obs_mod_rf<10 ),  # necessary to get useful bins with heatscatter()
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
  #   lines( c(-99,99), c(-99,99), col='red' )
  #   mtext( "FLUXCOM MTE-RF", line=1, adj=0.5 )

  #   df_dday_8d_agg <- df_dday_8d_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
    
  #   xvals <- fvarbins[1:nbins]+binwidth/2
  #   df_bins <- df_dday_8d_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_rf) ) %>% 
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


# #---------------------------------------------------------
# # MTE by site
# #---------------------------------------------------------
# for (sitename in unique(df_dday_mte_agg$mysitename)){
# 
#   sub <- filter( df_dday_mte_agg, mysitename == sitename )
# 
#   if (nrow(sub)>5){
#   
#     par( las=1, mar=c(4,4.5,2.5,0) )
#     with( 
#       filter( sub, ratio_obs_mod_mte<10 ),  # necessary to get useful bins with heatscatter()
#       plot( 
#         fvar, 
#         ratio_obs_mod_mte, 
#         xlab="fLUE",
#         ylab="GPP observed / GPP modelled",
#         xlim=c(0,1.2),
#         ylim=ylim,
#         cex=1.2,
#         pch=16,
#         col=add_alpha("black", 1),
#         main=""
#       ) 
#       
#     )
#     abline( h=1.0, lwd=0.5, lty=2 )
#     abline( v=1.0, lwd=0.5, lty=2 )
#     lines( c(-99,99), c(-99,99), col='red' )
#     mtext( sitename, line=1, adj=0.5 )
#     
#     sub <- sub %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
#     
#     xvals <- fvarbins[1:nbins]+binwidth/2
#     df_bins <- try( sub %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_mte) ) %>% 
#       summarise( peak=getpeak(ratio_obs_mod_mte), uhalfpeak=getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ) ) %>%
#       complete( infvarbin, fill = list( peak = NA ) ) %>% 
#       mutate( mids=xvals ))
#     
#     if (class(df_bins)!="try-error"){
#       rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
#       with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
#     }
#     
#   }
#   
# }



