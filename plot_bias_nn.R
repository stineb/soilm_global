.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

require(dplyr)
require(LSD)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )

getpeak <- function( vec ) {
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[ dens$y==max(dens$y) ]
  return( peak )   
}

getlhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x<peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
} 

getuhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x>peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
}


siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_modis_agg', 'df_dday_mte_agg', 

## load nice_agg to get data outside droughts
load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads nice_agg
load( "data/nice_nn_modis_agg_lue_obs_evi.Rdata" )   # loads modis_agg
load( "data/nice_nn_mte_agg_lue_obs_evi.Rdata" )   # loads mte_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )

# ##------------------------------------------------
# ## overwrite ratio_obs_mod again
# ##------------------------------------------------
# df_dday_agg$ratio_obs_mod       <- remove_outliers( df_dday_agg$ratio_obs_mod,       coef=5 )
# df_dday_modis_agg$ratio_obs_mod <- remove_outliers( df_dday_modis_agg$ratio_obs_mod, coef=5 )
# df_dday_mte_agg$ratio_obs_mod   <- remove_outliers( df_dday_mte_agg$ratio_obs_mod,   coef=5 )

##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------

## panel setup
magn <- 3
ncols <- 4
nrows <- 2
widths <- c(magn, 0.3*magn, magn, 0.3*magn )
heights <- 1.2*c(0.8*magn,0.8*magn,0.8*magn)
order <- matrix(seq(ncols*nrows),nrows,ncols,byrow=TRUE)

pdf( "fig/bias_vs_fvar_uncorrected.pdf", width=sum(widths), height=sum(heights) )

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
    par( las=1, mar=c(4,4.5,2.5,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,2.5)
    with( 
          filter( df_dday_agg, ratio_obs_mod_pmodel<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_pmodel, 
                      xlab="fLUE",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main=""
                    )
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "P-model", line=1, adj=0.5 )

    # ## add boxes for distribution within bins
    # nbins <- 10
    # binwidth <- 1.0/nbins
    # fvarbins <- seq( from=0, to=1, by=binwidth )
    # df_dday_agg <- df_dday_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
    # xvals <- fvarbins[1:nbins]+binwidth/2
    # df_bins <- df_dday_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_pmodel) ) %>% 
    #   summarise( peak=getpeak(ratio_obs_mod_pmodel), uhalfpeak=getuhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_pmodel, lev=0.75 ) ) %>%
    #   mutate( mids=xvals )
    # rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
    # with( df_bins, points( xvals, peak, pch='-', col="red", cex=5 ) )

    ## use only data during droughts for stats
    sub <- filter( df_dday_agg, is_drought_byvar==1 ) 
    stats <- analyse_modobs( sub$ratio_obs_mod_pmodel, sub$fvar, do.plot=FALSE )

    # # write stats into plot
    # x0 <- 0.05*xlim[2]
    # y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    # text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    # text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )
    
    ## draw the legend
    legend( "topleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", cex=0.8 )

    # Distribution 
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )
    # boxplot( filter( nice_agg, !is_drought_byvar )$ratio_obs_mod_pmodel, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50', xlim=c(0,5), at=0.5 )
    boxplot( filter( df_dday_agg, dday < 0 )$ratio_obs_mod_pmodel, outline=FALSE, axes=FALSE, ylim=ylim, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )

  #---------------------------------------------------------
  # MODIS
  #---------------------------------------------------------
    ## filtering is necessary to get sensible bins for coloring in the heatscatter function
    df_dday_modis_agg <- df_dday_modis_agg %>% filter( ratio_obs_mod_modis<5 )
    
    par( las=1, mar=c(4,4.5,2.5,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,3)
    with( 
          filter( df_dday_modis_agg, ratio_obs_mod_modis<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_modis, 
                      xlab="fLUE",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      cexplot=1.2,
                      main=""
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "MOD17A2H", line=1, adj=0.5 )

    ## add boxes for distribution within bins
    nbins <- 10
    binwidth <- 1.0/nbins
    fvarbins <- seq( from=0, to=1, by=binwidth )
    df_dday_modis_agg <- df_dday_modis_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_modis_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) ) %>% 
      summarise( peak=getpeak(ratio_obs_mod_modis), uhalfpeak=getuhalfpeak( ratio_obs_mod_modis, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_modis, lev=0.75 ),
                 q25=quantile(ratio_obs_mod_modis, probs=0.25), q75=quantile(ratio_obs_mod_modis, probs=0.75) ) %>%
      mutate( mids=xvals )
    
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
    
#     ## use only data during droughts for stats
#     sub <- filter( df_dday_modis_agg, is_drought_byvar==1 ) 
#     stats <- analyse_modobs( sub$ratio_obs_mod_modis, sub$fvar, do.plot=FALSE )
# 
#     ## write stats into plot
#     x0 <- 0.05*xlim[2]
#     y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
#     # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
#     # text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
#     text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( filter( df_dday_modis_agg, dday < 0 )$ratio_obs_mod_modis, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )

  #---------------------------------------------------------
  # MTE
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    with( 
          filter( df_dday_mte_agg, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with heatscatter()
          plot( 
                fvar, 
                ratio_obs_mod_mte, 
                xlab="fLUE",
                ylab="GPP observed / GPP modelled",
                xlim=c(0,1.2),
                ylim=ylim,
                cex=1.2,
                pch=16,
                col=add_alpha("black", 0.3),
                main=""
              ) 

        )
    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "FLUXCOM MTE", line=1, adj=0.5 )

    ## add boxes for distribution within bins
    nbins <- 10
    binwidth <- 1.0/nbins
    fvarbins <- seq( from=0, to=1, by=binwidth )
    df_dday_mte_agg <- df_dday_mte_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_mte_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_mte) ) %>% 
      summarise( peak=getpeak(ratio_obs_mod_mte), uhalfpeak=getuhalfpeak( ratio_obs_mod_mte, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_mte, lev=0.75 ),
      q25=quantile(ratio_obs_mod_mte, probs=0.25), q75=quantile(ratio_obs_mod_mte, probs=0.75) ) %>%
      mutate( mids=xvals )
    
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )
    

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( filter( df_dday_mte_agg, dday < 0 )$ratio_obs_mod_mte, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )


  #---------------------------------------------------------
  # MTE-RF
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    with( 
          filter( df_dday_mte_agg, ratio_obs_mod_rf<5 ),  # necessary to get useful bins with heatscatter()
          plot( 
                fvar, 
                ratio_obs_mod_rf, 
                xlab="fLUE",
                ylab="GPP observed / GPP modelled",
                xlim=c(0,1.2),
                ylim=ylim,
                cex=1.2,
                pch=16,
                col=add_alpha("black", 0.3),
                main=""
              ) 

        )
    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='red' )
    mtext( "FLUXCOM MTE-RF", line=1, adj=0.5 )

    nbins <- 10
    binwidth <- 1.0/nbins
    fvarbins <- seq( from=0, to=1, by=binwidth )
    df_dday_mte_agg <- df_dday_mte_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) )
    
    xvals <- fvarbins[1:nbins]+binwidth/2
    df_bins <- df_dday_mte_agg %>% group_by( infvarbin ) %>% filter( !is.na(infvarbin) & !is.na(ratio_obs_mod_rf) ) %>% 
      summarise( peak=getpeak(ratio_obs_mod_rf), uhalfpeak=getuhalfpeak( ratio_obs_mod_rf, lev=0.75 ), lhalfpeak=getlhalfpeak( ratio_obs_mod_rf, lev=0.75 ),
                 q25=quantile(ratio_obs_mod_rf, probs=0.25), q75=quantile(ratio_obs_mod_rf, probs=0.75) ) %>%
      mutate( mids=xvals )
    
    rect( xvals-0.02, df_bins$lhalfpeak, xvals+0.02, df_bins$uhalfpeak, col = add_alpha("white", 0.5) )
    with( df_bins, points( xvals, peak, pch='-', col="red", cex=2 ) )

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( filter( df_dday_mte_agg, dday < 0 )$ratio_obs_mod_rf, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )

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
#       filter( sub, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with heatscatter()
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
#     nbins <- 10
#     binwidth <- 1.0/nbins
#     fvarbins <- seq( from=0, to=1, by=binwidth )
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



