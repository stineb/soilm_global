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

plot_aligned_all <- function( ddf, ddf_8d, filn=NA ){
  ##------------------------------------------------
  ## Aligned plots for all sites combined into clusters (3 rows)
  ##------------------------------------------------
  require(dplyr)
  require(tidyr)
  source("../utilities/add_alpha.R")
  
  before <- 30
  after  <- 100

  if (!is.na(filn)) pdf( plotfiln, width=6, height=5 )

    ##--------------------------------------------------------
    ## fLUE
    ##--------------------------------------------------------
    par( las=1, mar=c(4,4,1,3), xpd=FALSE, xaxs="i", yaxs="r" )
    ylim <- c(0.0,1.5)
    xlim <- c(-20,80)
    plot( c(-before,after), ylim, type="n", xlab="Day after drought onset", ylab="fLUE / ratio", axes=FALSE, xlim=xlim )
    axis( 2, lwd = 1.75 )
    axis( 2, at = seq( ylim[1], ylim[2], by=0.1 ), labels = FALSE, tck=-0.01 )
    axis( 4, labels=FALSE, lwd = 1.75 )
    axis( 4, at = seq( ylim[1], ylim[2], by=0.1 ), labels = FALSE, tck=-0.01 )
    axis( 4 )
    axis( 1, xlab="days after drought onset", lwd=1.75 )
    axis( 1, at = seq( xlim[1], xlim[2], by=5 ), labels = FALSE, tck=-0.01 )
    abline( h=1.0, col='grey40', lwd=0.5 )
    
    axis(1,lwd=1.75);  axis(1,at=seq(xlim[1],xlim[2],by=20),labels=F,tck=-0.01)

    rect( 0, -99, after, 99, col=rgb(0,0,0,0.2), border=NA )

    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    median <- ddf %>% group_by( dday ) %>% 
                      summarise( fvar = median( fvar , na.rm=TRUE ) ) %>% 
                      complete( dday )
    upper  <- ddf %>% group_by( dday ) %>% 
                      summarise( fvar = quantile( fvar, 0.75, na.rm=TRUE ) ) %>% 
                      complete( dday )
    lower  <- ddf %>% group_by( dday ) %>% 
                      summarise( fvar = quantile( fvar, 0.25, na.rm=TRUE ) ) %>% 
                      complete( dday )

    polygon( c( median$dday, rev(median$dday) ), c( lower$fvar,  rev(upper$fvar) ),  col=add_alpha("black", 0.3), border=NA )
    lines( median, col='black', lwd=2 )


    ##--------------------------------------------------------
    ## Bias P-model
    ##--------------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    median <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = median( 1/bias_pmodel_norm , na.rm=TRUE ) ) %>% 
                      complete( dday )
    upper  <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = quantile( 1/bias_pmodel_norm, 0.75, na.rm=TRUE ) ) %>% 
                      complete( dday )
    lower  <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = quantile( 1/bias_pmodel_norm, 0.25, na.rm=TRUE ) ) %>% 
                      complete( dday )

    polygon( c( median$dday, rev(median$dday) ), c( lower$bias,  rev(upper$bias) ),  col=add_alpha("tomato", 0.3), border=NA )
    lines( median, col='tomato', lwd=2 )

    ##--------------------------------------------------------
    ## Bias BESS
    ##--------------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    median <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = median( 1/bias_bess_v1_norm , na.rm=TRUE ) ) %>% 
                      complete( dday )
    upper  <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = quantile( 1/bias_bess_v1_norm, 0.75, na.rm=TRUE ) ) %>% 
                      complete( dday )
    lower  <- ddf %>% group_by( dday ) %>% 
                      summarise( bias = quantile( 1/bias_bess_v1_norm, 0.25, na.rm=TRUE ) ) %>% 
                      complete( dday )

    polygon( c( median$dday, rev(median$dday) ), c( lower$bias,  rev(upper$bias) ),  col=add_alpha("royalblue3", 0.3), border=NA )
    lines( median, col='royalblue3', lwd=2 )

    ##--------------------------------------------------------
    ## Bias MODIS
    ##--------------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    median <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = median( 1/bias_modis_norm , na.rm=TRUE ) ) %>% 
                         complete( dday )
    upper  <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = quantile( 1/bias_modis_norm, 0.75, na.rm=TRUE ) ) %>% 
                         complete( dday )
    lower  <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = quantile( 1/bias_modis_norm, 0.25, na.rm=TRUE ) ) %>% 
                         complete( dday )

    polygon( c( median$dday, rev(median$dday) ), c( lower$bias,  rev(upper$bias) ),  col=add_alpha("springgreen1", 0.3), border=NA )
    lines( median, col='springgreen1', lwd=2 )

    ##--------------------------------------------------------
    ## Bias VPM
    ##--------------------------------------------------------
    ## Get median level of the three variables within each bin, pooling data for all days and instances (drought events)
    median <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = median( 1/bias_vpm_norm , na.rm=TRUE ) ) %>% 
                         complete( dday )
    upper  <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = quantile( 1/bias_vpm_norm, 0.75, na.rm=TRUE ) ) %>% 
                         complete( dday )
    lower  <- ddf_8d %>% group_by( dday ) %>% 
                         summarise( bias = quantile( 1/bias_vpm_norm, 0.25, na.rm=TRUE ) ) %>% 
                         complete( dday )

    polygon( c( median$dday, rev(median$dday) ), c( lower$bias,  rev(upper$bias) ),  col=add_alpha("springgreen4", 0.3), border=NA )
    lines( median, col='springgreen4', lwd=2 )

    legend( "bottomleft", c("fLUE", "P-model", "BESS", "MODIS", "VPM"), bty="n", lty=1, lwd=2, col=c( "black" ,"tomato", "royalblue3","springgreen1", "springgreen4"), cex=1.0 )

  if (!is.na(filn)) dev.off()

}

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
load( "data/data_aligned_agg.Rdata" )
plot_aligned_all( df_dday_agg, df_dday_8d_agg, filn=NA )

