lmp <- function(modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

inverse <- function( x ){ 1/x }

normalise <- function( x, norm ){ x/norm }

smooth_runminmax <- function( x, y ){  
  source( "../utilities/cutna_headtail.R" )
  idxs_drop <- cutna_headtail( y )
  if (length(idxs_drop)>0) { x <- x[-idxs_drop] }
  if (length(idxs_drop)>0) { y <- y[-idxs_drop] }
  x[is.nan(x)] <- NA
  x[is.infinite(x)] <- NA
  fld <- approx( x, y, xout=x )$y
  min <- smooth.spline( x, caTools::runmin( fld, k=7 ), spar=0.01 )$y
  max <- smooth.spline( x, caTools::runmax( fld, k=7 ), spar=0.01 )$y
  tmp <- approx( x, zoo::rollmean( fld, k=7, fill=NA ), xout=x )$y
  idxs<- which(!is.na(tmp))
  tmp <- smooth.spline( x[idxs], tmp[idxs], spar=0.01 )$y
  mean <- rep(NA, length(x))
  mean[idxs] <- tmp
  return( list( x=x, min=min, max=max, mean=mean ) )
}

spline_with_gaps <- function( xvals, yvals, nice, spar=NULL ){
  idxs <- which( !is.na(yvals) & !is.nan(yvals)  )
  tmp <- smooth.spline( xvals[idxs], yvals[idxs], spar=spar )
  tmp <- data.frame( year_dec=tmp$x, yvals=tmp$y )
  nice <- nice %>% left_join( tmp, by="year_dec" )
  return( nice$yvals )       
}


plot_bysite_tseries <- function( sitename, makepdf=TRUE ){

  require(dplyr)
  require(lubridate)

  # ##--------------------
  # ## Load stuff
  # ##--------------------
  # sitename = "AR-Vir"
  # makepdf = TRUE

  ## filter for this site
  nice    <- filter( nice_agg, mysitename==sitename )
  nice_8d <- filter( nice_8d_agg, mysitename==sitename )
  df_dday_aggbydday <- filter( df_dday_aggbydday_agg, mysitename==sitename )
  df_dday_aggbydday_8d <- filter( df_dday_aggbydday_8d_agg, mysitename==sitename )

  ## get additional data
  infil <- paste0( myhome, "/data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_lue_obs_evi.Rdata" ) 
  load( infil ) ## gets list 'nn_fluxnet'
  minmax           <- nn_fluxnet[[ sitename ]]$minmax          
  droughts         <- nn_fluxnet[[ sitename ]]$droughts        
  ##--------------------

  # require( dplyr )
  # require( tidyr )
  # require( minpack.lm )
  # require( LSD )
  # require( cluster )
  # require( broom )
  # require( zoo )
  require(stats)
  require(readr)
  require(lubridate)

  source("../utilities/get_consecutive.R")
  source( "../utilities/add_alpha.R" )

  siteinfo <- read_csv( "siteinfo_fluxnet2015_sofun.csv" )

  fapar_data <- "fpar"
  use_fapar <- FALSE
  char_fapar <- "_withFPAR"

  nice    <- nice %>% mutate( year_dec = decimal_date(date) )
  nice_8d <- nice_8d %>% mutate( year_dec = decimal_date(date) )

  # ##------------------------------------------------
  # ## PLOT TIME SERIES
  # ##------------------------------------------------
  #   filn <- paste0("fig/plot_bysite_", sitename, "_tseries.pdf")  

  #   nyears <- ceiling(range(nice$year_dec)[2]) - floor(range(nice$year_dec)[1])

  #   magn <- 3.5
  #   ncols <- 1

  #   ## with 2 panels only
  #   heights <- magn
  #   nrows <- 1

  #   widths <- 10
  #   # widths  <- rep( nyears / 5,ncols) * magn

  #   if (makepdf) pdf( filn, width=sum(widths), height=sum(heights), bg="white" )

  #     panel <- layout(
  #                     matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
  #                     # matrix( order, nrows, ncols, byrow=TRUE ),
  #                     widths=widths,
  #                     heights=heights,
  #                     TRUE
  #                     )
  #     # layout.show( panel )

  #     cols <- c( "black", "royalblue3", "tomato", "springgreen3", "orchid" )
  #     varnam <- "fLUE"
  #     descr_outputset <- c(
  #       "NN, good days (LUE ~ temperature + VPD + PPFD)",
  #       "NN, all days (LUE ~ temperature + VPD + PPFD + soil moisture)"
  #       )
  #     lue <- TRUE
  #     data_fapar <- "fpar"

  #     xlim <- range(nice$date)

  #     # ## Legend
  #     # legend( "topright", c("observed", expression( paste( "NN"[pot])), expression(paste("NN"[act])) ), bty="n", lty=1, lwd=2, col=c("grey35", cols[1], cols[2]), cex=0.9 )

  #     ##----------------------------------------------------------------------------------------
  #     ## 2nd panel: fLUE and EVI
  #     ##----------------------------------------------------------------------------------------
  #     par( las=1, mar=c( 4, 4.4, 1, 3 ), xpd=FALSE )
  #     with( nice, plot( year_dec, fvar, type="n", xlab="year", ylab="unitless", axes=FALSE, ylim=c(0,1.6) ) )
  #     axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )),     labels=FALSE )
  #     axis( 1, at=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) ))+0.5, labels=unique( c( floor(nice$year_dec)[1], ceiling(nice$year_dec) )), tck=0.0 )
  #     axis( 2 )
  #     abline( h=1, lwd=0.5 )

  #     ## rectangles for droughts
  #     if (!is.null(droughts) && nrow(droughts)>0 ){
  #       rect( nice$year_dec[droughts$idx_start], rep( -99, nrow(droughts) ), nice$year_dec[(droughts$idx_start+droughts$len-1)], rep( 99, nrow(droughts) ), col=rgb(0,0,0,0.2), border=NA )
  #     }

  #     ## solid line for fLUE (smoothed)
  #     with( nice, lines( year_dec, fvar_smooth, col=cols[1], lwd=1 ) )

  #     ## Uncertainty range of fLUE
  #     polygon( c( minmax$year_dec, rev(minmax$year_dec) ), c( smooth.spline( minmax$year_dec, minmax$fvar_min_filled, spar=0.01 )$y, rev( smooth.spline( minmax$year_dec, minmax$fvar_max_filled, spar=0.01 )$y ) ), border=NA, col=add_alpha(cols[1],0.5) )

  #     ## Ratio P-model
  #     with( nice, lines( year_dec, spline_with_gaps( year_dec, ratio_obs_mod_pmodel, nice, spar=0.01 ), col=cols[3], lwd=1.5 ))

  #     ## Ratio BESS
  #     with( nice, lines( year_dec, spline_with_gaps( year_dec, ratio_obs_mod_bess_v1, nice, spar=0.01 ), col=cols[2], lwd=1.5 ))

  #     ## Ratio MODIS
  #     with( nice_8d, lines( year_dec, spline_with_gaps( year_dec, ratio_obs_mod_modis, nice_8d, spar=0.08 ), col=cols[4], lwd=1.5 ))

  #     ## Ratio VPM
  #     with( nice_8d, lines( year_dec, spline_with_gaps( year_dec, ratio_obs_mod_vpm, nice_8d, spar=0.08 ), col=cols[5], lwd=1.5 ))
      
  #     # ## Plot time series of EVI
  #     # lines( nice$year_dec, nice[[ fapar_data ]], col="springgreen3", lwd=2 )
  #     # # if (!missing_modis_fpar){
  #     # #   lines( modis_fpar$year_dec, modis_fpar$fapar_modis*1e1, col='springgreen3', lwd=1 )
  #     # # }
  #     # 
  #     # # ## Add rectangle for fAPAR extremes in time series
  #     # # if (!is.null(fapar_extremes)){
  #     # #   # abline( v=fapar_extremes$year_dec_start, col='red' )
  #     # #   # abline( v=fapar_extremes$year_dec_end, col='red' )
  #     # #   for (idx in 1:nrow(fapar_extremes)){
  #     # #     idx_start <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_start[idx]) )
  #     # #     idx_end   <- which.min( abs(nice$year_dec-fapar_extremes$year_dec_end[idx]) )
  #     # #     if (idx_start!=1 && idx_end!=1 ){
  #     # #         rect( nice$year_dec[ idx_start ], ylim[1], nice$year_dec[ idx_end ], ylim[2], col=add_alpha('darkgoldenrod', 0.3), border=NA, lwd=2 )
  #     # #     }
  #     # #   }
  #     # # }
  #     # 
  #     # ## Legend
  #     # # legend( "bottomright", c("fLUE", "MODIS EVI","MODIS FPAR"), bty="n", lty=1, lwd=c(1,2,1), col=c("tomato", "palegreen3", "springgreen3"), cex=1.2 )
  #     # legend( "bottomright", c("fLUE", "MODIS EVI"), bty="n", lty=1, lwd=2, col=c("tomato", "springgreen3"), cex=0.9 )

  #   if (makepdf) dev.off()


  ##------------------------------------------------
  ## PLOT ALIGNED
  ##------------------------------------------------
  filn <- paste0("fig/plot_bysite_", sitename, "_aligned.pdf")  

  before <- 30
  after  <- 100
  alpha <- 0.3/(nrow(droughts))
  xvals <- (-before:after)+1
 
  magn <- 2.4
  ncols <- 1
  nrows <- 1
  heights <- 1.4*magn
  widths  <- rep(1.8,ncols)*magn

  if (makepdf) pdf( filn, width=sum(widths), height=sum(heights), bg="white" )

    panel <- layout(
                    matrix(c(1:(nrows*ncols)),nrows,ncols,byrow=TRUE),
                    # matrix( order, nrows, ncols, byrow=TRUE ),
                    widths=widths,
                    heights=heights,
                    TRUE
                    )
    # layout.show( panel )

    ##--------------------------------------------------------
    ## fLUE
    ##--------------------------------------------------------
      xvals <- (-before:after)+1

      par( las=1, mar=c(4,4,1,4), xpd=FALSE )
      plot( c(-before,after), c(0,1.6), type="n", xlab="day after drought onset", ylab="soil water content (fraction)", axes=FALSE ) # , col.lab="blue"
      axis( 1, xlab="days after drought onset" )
      axis( 2 )
      abline( h=1.0, col='grey40', lwd=0.5 )

      rect( 0, -99, droughts$len, 99, col=rgb(0,0,0,alpha), border=NA )

      ## plot polygon for fvar
      upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_upp, xout=xvals )$y
      lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_low, xout=xvals )$y
      mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$fvar_med, xout=xvals )$y
      idxs  <- which( !is.na(upper) & !is.na(lower) )

      polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha("tomato", 0.4), border=NA )
      lines( xvals, mid, col="tomato" )


    ##--------------------------------------------------------
    ## P-model bias 
    ##--------------------------------------------------------
      par( new=TRUE )
      xvals <- (-before:after)+1

      if ("bias_pmodel_med" %in% names(df_dday_aggbydday) && any(!is.na(df_dday_aggbydday$bias_pmodel_upp))){
        
        # norm <- mean( 1/df_dday_aggbydday$bias_pmodel_med[1:30], na.rm=TRUE )
        norm <- 1
        
        upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_pmodel_upp, xout=xvals )$y %>% inverse() %>% normalise( norm )
        lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_pmodel_low, xout=xvals )$y %>% inverse() %>% normalise( norm )
        mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_pmodel_med, xout=xvals )$y %>% inverse() %>% normalise( norm )
        idxs  <- which( !is.na(upper) & !is.na(lower) )

        ## plot polygon
        polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('orchid', 0.3), border=NA )
        lines( xvals, mid, col=add_alpha('orchid', 1), lwd=1 )

      }

    ##--------------------------------------------------------
    ## BESS bias 
    ##--------------------------------------------------------
      par( new=TRUE )
      xvals <- (-before:after)+1

      if ("bias_bess_v1_med" %in% names(df_dday_aggbydday) && any(!is.na(df_dday_aggbydday$bias_bess_v1_upp))){

        # norm <- mean( 1/df_dday_aggbydday$bias_bess_v1_med[1:30], na.rm=TRUE )
        norm <- 1
        
        upper <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_bess_v1_upp, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        lower <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_bess_v1_low, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        mid   <- approx( df_dday_aggbydday$dday, df_dday_aggbydday$bias_bess_v1_med, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        idxs  <- which( !is.na(upper) & !is.na(lower) )

        ## plot polygon
        polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('royalblue3', 0.3), border=NA )
        lines( xvals, mid, col=add_alpha('royalblue3', 1), lwd=1 )

      }

    ##--------------------------------------------------------
    ## MODIS bias 
    ##--------------------------------------------------------
      par( new=TRUE )
      xvals <- (-before:after)+1

      # norm <- mean( 1/df_dday_aggbydday_8d$bias_modis_med[1:4], na.rm=TRUE )
      norm <- 1
      
      upper <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_modis_q75, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
      lower <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_modis_q25, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
      mid   <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_modis_med, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
      idxs  <- which( !is.na(upper) & !is.na(lower) )

      ## plot polygon
      polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('springgreen3', 0.3), border=NA )
      lines( xvals, mid, col=add_alpha('springgreen3', 1), lwd=1 )

    ##--------------------------------------------------------
    ## VPM bias 
    ##--------------------------------------------------------
      par( new=TRUE )
      xvals <- (-before:after)+1

      if ("bias_vpm_med" %in% names(df_dday_aggbydday_8d) && any(!is.na(df_dday_aggbydday_8d$bias_vpm_med))){

        # norm <- mean( 1/df_dday_aggbydday_8d$bias_vpm_med[1:4], na.rm=TRUE )
        norm <- 1
        
        upper <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_vpm_q75, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        lower <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_vpm_q25, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        mid   <- approx( df_dday_aggbydday_8d$dday*8, df_dday_aggbydday_8d$bias_vpm_med, xout=xvals )$y %>% inverse() %>% normalise( norm ) 
        idxs  <- which( !is.na(upper) & !is.na(lower) )

        ## plot polygon
        polygon( c( xvals[idxs], rev(xvals[idxs])), c( lower[idxs], rev(upper[idxs])), col=add_alpha('springgreen4', 0.3), border=NA )
        lines( xvals, mid, col=add_alpha('springgreen4', 1), lwd=1 )

      }

      # text( -before, 1.2, sitename, cex=1.2, adj=c(0,0), font=2 )
      # text( -before, 1.1, dplyr::filter( siteinfo, siteinfo$mysitename==sitename )$classid, cex=1.2, adj=c(0,0), font=1 )

      # legend( "bottomleft", c("soil water content (fraction)", expression( paste( "VPD"^{-1}, " (relative change)" ) ) ), bty="n", lty=1, lwd=2, col=c("blue", "steelblue3"), cex=1.0 )


  if (makepdf) dev.off()


}

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 )$mysitename

load( "data/nice_nn_agg_lue_obs_evi.Rdata" )
load( "data/nice_nn_8d_agg_lue_obs_evi.Rdata" )
load( "data/data_aligned_agg.Rdata" )

do.sites <- "FR-Pue"

print("creating time series plots for all sites ...")
for (sitename in do.sites){
  plot_bysite_tseries( sitename, makepdf=TRUE )
}
print("... done.")

