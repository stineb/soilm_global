syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

require(dplyr)
require(LSD)

source( "analyse_modobs.R" )
source( "remove_outliers.R" )
source( "compl_df_flue_est.R" )

getpeak <- function( vec ) {
  vec  <- vec[!is.na(vec)] %>% remove_outliers( coef=5.0 )
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[ dens$y==max(dens$y) ]
  return( peak )   
}

getlhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  vec  <- vec[!is.na(vec)] %>% remove_outliers( coef=5.0 )
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x<peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
} 

getuhalfpeak <- function( vec, lev=0.5 ) {
  require(dplyr)
  vec  <- vec[!is.na(vec)] %>% remove_outliers( coef=5.0 )
  dens <- density( vec, kernel=c("gaussian") )
  peak <- dens$x[dens$y==max(dens$y)]
  df_tmp <- tibble( x=dens$x, y=dens$y ) %>% filter( x>peak )
  halfpeak <- df_tmp$x[ which.min( abs( df_tmp$y - lev * max(dens$y) ) ) ]
  return( halfpeak )   
}


siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aligned aggregated data
load( "data/data_aligned_agg.Rdata" ) # loads 'df_dday_agg', 'df_dday_8d_agg', 'df_dday_mte_agg', 'df_dday_bess_agg', 'df_dday_vpm_agg'

# ## Estimate soil moisture correction (adds column 'flue_est' to dataframe)
# load( "data/linearfit2_ratio.Rdata" )
# df_dday_8d_agg    <- compl_df_flue_est( df_dday_8d_agg, linearfit2, x0_fix=0.9  )

## group data by fLUE
nbins <- 10
binwidth <- 1.0/nbins
fvarbins <- seq( from=0, to=1, by=binwidth )
xvals <- fvarbins[1:nbins]+binwidth/2

df_dday_8d_agg    <- df_dday_8d_agg %>% mutate( infvarbin = cut( fvar, breaks = fvarbins ) ) %>%
                                        mutate( ifelse( is.nan(ratio_obs_mod_pmodel), NA, ratio_obs_mod_pmodel ) )

##------------------------------------------------
## GPPobs/GPPmod vs. fLUE
##------------------------------------------------

## panel setup
magn <- 3
ncols <- 2
nrows <- 3
widths <- c(magn, 0.92*magn )
heights <- 1.2*c(0.6*magn,0.6*magn,0.7*magn)
order <- matrix(seq(ncols*nrows),nrows,ncols,byrow=TRUE)

pdf( "fig/bias_vs_fvar_uncorrected_8d.pdf", width=sum(widths), height=sum(heights) )

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
    par( las=1, mar=c(2,4.5,2.5,1) )
    xlim <- c(0,1.2)
    ylim <- c(0,2.5)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_pmodel<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_pmodel, 
                      xlab="",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main="",
                      cexplot = 1.2
                    )
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "P-model", line=0.5, adj=0, font=2, cex=0.8 )
    
    ## draw the legend
    legend( "topleft", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray60", "navy", "red", "yellow"))(5), bty="n", cex=0.8 )
    
  #---------------------------------------------------------
  # MODIS
  #---------------------------------------------------------
    par( las=1, mar=c(2,2.5,2.5,1) )
    xlim <- c(0,1.2)
    ylim <- c(0,3)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_modis<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
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
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "MOD17A2H", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # BESS v1
  #---------------------------------------------------------
    par( las=1, mar=c(2,4.5,2.5,1) )
    xlim <- c(0,1.2)
    ylim <- c(0,3)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_bess_v1<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_bess_v1, 
                      xlab="",
                      ylab="GPP observed / GPP modelled",
                      xlim=xlim,
                      ylim=ylim,
                      main="",
                      cexplot = 1.2
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "BESS v1", line=0.5, adj=0, font=2, cex=0.8 )
    

  #---------------------------------------------------------
  # BESS v2
  #---------------------------------------------------------
    par( las=1, mar=c(2,2.5,2.5,1) )
    xlim <- c(0,1.2)
    ylim <- c(0,3)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_bess_v2<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_bess_v2, 
                      xlab="",
                      ylab="",
                      xlim=xlim,
                      ylim=ylim,
                      main="",
                      cexplot = 1.2
                    ) 
        )

    abline( h=1.0, lwd=0.5, lty=2 )
    abline( v=1.0, lwd=0.5, lty=2 )
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "BESS v2", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # VPM
  #---------------------------------------------------------
    par( las=1, mar=c(4,4.5,2.5,1) )
    xlim <- c(0,1.2)
    ylim <- c(0,3)
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_vpm<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_vpm, 
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
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "VPM", line=0.5, adj=0, font=2, cex=0.8 )

  #---------------------------------------------------------
  # MTE
  #---------------------------------------------------------
    par( las=1, mar=c(4,2.5,2.5,1) )
    with( 
          filter( df_dday_8d_agg, ratio_obs_mod_mte<5 ),  # necessary to get useful bins with heatscatter()
          heatscatter( 
                      fvar, 
                      ratio_obs_mod_mte, 
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
    lines( c(-99,99), c(-99,99), col='black' )
    mtext( "FLUXCOM MTE", line=0.5, adj=0, font=2, cex=0.8 )


dev.off()


