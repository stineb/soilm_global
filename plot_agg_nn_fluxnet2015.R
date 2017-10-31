.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

library(dplyr)
library(LSD)

source( paste( myhome, "sofun/utils_sofun/analysis_sofun/analyse_modobs.R", sep="" ) )
source( paste( myhome, "sofun/utils_sofun/analysis_sofun/remove_outliers.R", sep="" ) )

## Manual settings ----------------
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
nrep       = 5
dotrain    = FALSE
overwrite_modis = TRUE
overwrite_mte = TRUE
##--------------------------

## check and override if necessary
if ( nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
  plotlue <- TRUE
  if (nam_target=="lue_obs_evi"){
    fapar_data <- "evi"
  } else if (nam_target=="lue_obs_fpar"){
    fapar_data <- "fpar"
  }
  if (use_fapar){
    print("WARNING: setting use_fapar to FALSE")
    use_fapar <- FALSE
  }
}

## identifier for output files
if (use_fapar){
  if (nam_target=="lue_obs_evi"){
    char_fapar <- "_withEVI"
  } else if (nam_target=="lue_obs_fpar"){
    char_fapar <- "_withFPAR"
  } else {
    print("ERROR: PROVIDE VALID FAPAR DATA!")
  }
} else {
  char_fapar <- ""
}

if (use_weights){
  char_wgt <- "_wgt"
} else {
  char_wgt <- ""
}

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( paste( "data/nice_agg_", char_wgt, nam_target, ".Rdata", sep="" ) )       # loads 'nice_agg'
load( paste( "data/nice_resh_", char_wgt, nam_target, ".Rdata", sep="" ) )      # loads 'nice_resh'
load( paste( "data/nice_modis_agg_", char_wgt, nam_target, ".Rdata", sep="" ) ) # loads 'nice_to_modis_agg'
load( paste( "data/nice_mte_agg_", char_wgt, nam_target, ".Rdata", sep="" ) )   # loads 'nice_to_mte_agg'

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_modis_agg', 'df_dday_mte_agg', 

## add vegetation type info to nice_agg
nice_agg <- nice_agg %>% left_join( dplyr::select( siteinfo, mysitename, classid ), by="mysitename" )


##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------
  df_dday_agg$ratio_obs_mod             <- remove_outliers( df_dday_agg$ratio_obs_mod, coef=5 )
  df_dday_modis_agg$ratio_obs_mod_modis <- remove_outliers( df_dday_modis_agg$ratio_obs_mod_modis, coef=5 )
  df_dday_mte_agg$ratio_obs_mod_mte     <- remove_outliers( df_dday_mte_agg$ratio_obs_mod_mte, coef=5 )

  magn <- 3
  ncols <- 4
  nrows <- 2
  widths <- c(magn, 0.3*magn, magn, 0.3*magn )
  heights <- 1.2*c(magn,0.8*magn)
  order <- matrix(c(1,1,2,2,3,4,5,6),nrows,ncols,byrow=TRUE)

  pdf( "fig_nn_fluxnet2015/bias_vs_fvar/bias_vs_fvar_ALL.pdf", width=sum(widths), height=sum(heights) )

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
    par( las=1, mar=c(5,4.5,3,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,2.5)
    with( 
          dplyr::filter( df_dday_agg, ratio_obs_mod<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod, 
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

    ## use only data during droughts for stats
    sub <- dplyr::filter( df_dday_agg, is_drought_byvar==1 ) 
    stats <- analyse_modobs( sub$ratio_obs_mod, sub$fvar, do.plot=FALSE )

    # write stats into plot
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )
    
    ## draw the legend
    legend( "topleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", inset=c(0.05,0.25), cex=0.8 )

    # Distribution 
    par( las=1, mar=c(5,0,3,5), xpd=FALSE )

    boxplot( dplyr::filter( nice_agg, !is_drought_byvar )$ratio_obs_mod, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50', xlim=c(0,5), at=0.5 )
    lines( c(-2,1), c(1,1), lwd=0.5, lty=2 )
    # abline( h=1.0, lwd=0.5, lty=2 )


    #---------------------------------------------------------
    # MODIS
    #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    xlim <- c(0,1.2)
    ylim <- c(0,4)
    with( 
          dplyr::filter( df_dday_modis_agg, ratio_obs_mod_modis<5 ),  # necessary to get useful bins with heatscatter()
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
    mtext( "MODIS", line=1, adj=0.5 )

    ## use only data during droughts for stats
    sub <- dplyr::filter( df_dday_modis_agg, is_drought_byvar==1 ) 
    stats <- analyse_modobs( sub$ratio_obs_mod_modis, sub$fvar, do.plot=FALSE )

    ## write stats into plot
    x0 <- 0.05*xlim[2]
    y0 <- 0.8*(ylim[2]-ylim[1])+ylim[1]      
    # text( x0, y0, paste( "RMSE =", format( stats$rmse, digits = 2 ), " (", format( stats$prmse, digits = 2 ), "%)", sep="" ), adj=0.0, cex=0.8 )
    # text( x0, y0+0.15, bquote( italic(R)^2 == .(format( stats$rsq, digits = 2) ) ),  adj=0.0, cex=0.8 )
    text( x0, y0+0.3, paste( "N =", format( stats$N, digits = 1 ) ), adj=0.0, cex=0.8 )

    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( 1.0 / filter( nice_to_modis_agg, !is_drought_byvar )$bias_modis, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )


    #---------------------------------------------------------
    # MTE
    #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,0) )
    with( 
          dplyr::filter( df_dday_mte_agg, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with heatscatter()
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


    #---------------------------------------------------------
    # Distribution 
    #---------------------------------------------------------
    par( las=1, mar=c(4,0,2.5,2), xpd=FALSE )

    boxplot( 1.0 / filter( nice_to_mte_agg, !is_drought_byvar )$bias_mte, outline=FALSE, ylim=ylim, axes=FALSE, col='grey50' )
    abline( h=1.0, lwd=0.5, lty=2 )

  dev.off()





