require( fields, quietly = TRUE )
require( sp, quietly = TRUE )
require( maptools, quietly = TRUE )
require( dplyr, quietly = TRUE )  

source("../utilities/add_alpha.R")
source("../utilities/mycolorbar.R")

## load global GPP time series, prepared by plot_effects_gpp_tseries.R
load("data/gpp_glob_tseries.Rdata")

## load fields, derived by map_effects_gpp_mean.R
load("data/gpp_loss.Rdata" )

##------------------------------------------------------------------------
## Panel
##------------------------------------------------------------------------
magn <- 4
ncols <- 3
nrows <- 2
widths <- rep(1.4*magn,ncols)
widths[2] <- 0.15*widths[1]
widths[3] <- 0.80*widths[1]
heights <- rep(magn,nrows)
order <- matrix( c(1:6), nrows, ncols, byrow=TRUE )

pdf( "fig/fig_2.pdf", width=sum(widths), height = sum(heights) )
  
  par(las=1)

  panel <- layout(
          order,
          widths=widths,
          heights=heights,
          TRUE
          )
  # layout.show( panel )

  ##------------------------------------------------------------------------
  ## map GPP loss
  ##------------------------------------------------------------------------
  # out <- gpp_s0 - gpp_s1b
  out <- (1-(gpp_s1b / gpp_s0))*100
  lon <- seq(-179.75, 179.75, 0.5)
  lat <- seq(-89.75, 89.75, 0.5)
  
  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)
  
  a <- sapply( lat.labels, function(x) if (x>0) {bquote(.(x)*degree ~N)} else if (x==0) {bquote(.(x)*degree)} else {bquote(.(-x)*degree ~S)} )
  b <- sapply( lon.labels, function(x) if (x>0) {bquote(.(x)*degree ~E)} else if (x==0) {bquote(.(x)*degree)} else {bquote(.(-x)*degree ~W)})
  
  color <- c( "wheat", "tomato2", "tomato4" )
  
  # lev = c( 0, 400, 10 )
  lev=c( 0, 70, 7 )
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=max(out, na.rm=TRUE), minval=min(out, na.rm=TRUE) )
  
  par( mar=c(4,3,2,1),xaxs="i", yaxs="i",las=1)
  image(
          lon, lat,
          out,
          ylim=c(-60,85),
          yaxt="n", xaxt="n",
          col=out.mycolorbar$colors, breaks=out.mycolorbar$margins,
          xlab="", ylab=""
          )
  map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
  
  axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.9, lwd=1.5 )
  axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 4, at=lat.labels, lab=F, lwd=1.5 )
  axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.9, lwd=1.5 )
  axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 3, at=lon.labels, lab=F, lwd=1.5 )
  axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )
  
  mtext( "a)", font=2, adj = 0, line = 0.5, cex = 1.2 )

  ## Color key
  par( mar=c(4,3,2,1),xaxs="i", yaxs="i",las=1)
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1900 )

  ##------------------------------------------------------------------------
  ## time series
  ##------------------------------------------------------------------------
  par( mar=c(4,4.2,2,1), mgp =c(2.8,1,0))
  with( df, plot( year, gpp_s1b, type="n", ylim=c(100,160), xlab="Year", ylab=expression( paste("Global GPP (Pg C yr"^-1, ")" ) ), cex.lab=1.2 ) )
  with( df,  polygon( c( year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), border = NA, col=rgb(0,0,0,0.3) ) )

  with( df, lines( year, gpp_s1b, lty=2, lwd=1.5, col="tomato" ) )
  with( df, lines( year, gpp_s0, col="tomato", lwd=1.5 ) )
  with( df, lines( year, gpp_modis, col="orchid", lwd=1.5 ) )
  with( df, lines( year, gpp_vpm, col="springgreen3", lwd=1.5 ) )
  with( df, lines( year, gpp_bess, col="royalblue3", lwd=1.5 ) )
  with( df, lines( year, gpp_mte, col="darkgoldenrod3", lwd=1.5 ) )

  legend("topleft", c( "P-model, s0", "P-model, s1b (grey range: s1a - s1c)", "MODIS", "VPM", "BESS", "MTE" ), bty = "n", lty = c(1,2,1,1,1,1), lwd=1.5, col=c("tomato", "tomato", "orchid", "springgreen3", "royalblue3", "darkgoldenrod3") )
  mtext( "b)", font=2, adj = 0, line = 0.5, cex = 1.2 )

  ##------------------------------------------------------------------------
  ## map GPP loss trend
  ##------------------------------------------------------------------------
  # plot_map( trend_rel[1,,], lev = c(-0.5, 0.5, 10), positive = FALSE, minval = -1.4, maxval = 2.0 )
  out <- trend_rel[1,,]
  out[which( trend_rel[2,,] == 0 )] <- NA
  color <- rev(c( "royalblue4", "royalblue2", "wheat", "tomato2", "tomato4" ))

  lev = c(-0.5, 0.5, 10)
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, minval = -1.4, maxval = 2.0 )
  
  par( mar=c(4,3,2,1),xaxs="i", yaxs="i",las=1)
  image(
          lon, lat,
          out,
          ylim=c(-60,85),
          yaxt="n", xaxt="n",
          col=out.mycolorbar$colors, breaks=out.mycolorbar$margins,
          xlab="", ylab=""
          )
  map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
  
  axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.9, lwd=1.5 )
  axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 4, at=lat.labels, lab=F, lwd=1.5 )
  axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.9, lwd=1.5 )
  axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 3, at=lon.labels, lab=F, lwd=1.5 )
  axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )

  # ## add stippling
  # incl <- which( trend_rel[2,,] == 1 )
  # grd <- expand.grid( x=lon, y=lat )
  # incl <- incl[ seq(4,length(incl), by=4) ]
  # points( grd$x[incl], grd$y[incl], pch=".", cex=1 )
  # 
  mtext( "c)", font=2, adj = 0, line = 0.5, cex = 1.2 )

  ## Color key
  par( mar=c(4,3,2,1), xaxs="i", yaxs="i", las=1)
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, minval = -1.4, maxval = 2.0 )


  ##------------------------------------------------------------------------
  ## time series GPP loss %
  ##------------------------------------------------------------------------
  linmod_s1b_releff <- lm( -100*(1-gpp_s1b/gpp_s0) ~ year, data=df )

  ci <- confint(linmod_s1b_releff)
  newx <- seq( min(df$year)-5, max(df$year)+5, length.out=100 )
  preds <- predict( linmod_s1b_releff, newdata = data.frame( year = newx ), interval = 'confidence' )

  par( mar=c(4,4.2,2,1), mgp =c(2.8,1,0))
  with( df, plot( year, -100*(1-gpp_s1b/gpp_s0), pch=16, col="black", lty=1, xlab="Year", ylim=c(-16, -14), ylab=expression( paste("Reduction in global GPP (%)" ) ), xlim=c(1981,2017), cex.lab=1.2 ) )
  polygon( c( rev(newx), newx ), c( rev(preds[,3]), preds[,2] ), col = rgb(0,0,0,0.3), border = NA )
  abline( linmod_s1b_releff, lwd = 1.5, col="red" )
  text( 1984, -15.9, bquote( italic("slope") ==  .(format(coef(linmod_s1b_releff)[2], digits = 2) ) ~ "["* .(format( ci[2,1], digits = 2 )) ~ "-" ~ .(format( ci[2,2], digits = 2 ))*"]" ~ "yr" ^-1 ), adj=0, cex=1.2 )
  mtext( "d)", font=2, adj = 0, line = 0.5, cex = 1.2 )
  
dev.off()