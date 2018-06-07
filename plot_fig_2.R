require( fields, quietly = TRUE )
require( sp, quietly = TRUE )
require( maptools, quietly = TRUE )
require( dplyr, quietly = TRUE )  

source("../utilities/add_alpha.R")
source("../utilities/mycolorbar.R")

## load global GPP time series, prepared by plot_effects_gpp_tseries.R
load("data/gpp_glob_tseries.Rdata")

## load fields
load("data/gpp_loss.Rdata" )

##------------------------------------------------------------------------
## Panel
##------------------------------------------------------------------------
magn <- 4
ncols <- 3
nrows <- 1
widths <- rep(1.4*magn,ncols)
widths[2] <- 0.17*widths[1]
widths[3] <- 0.80*widths[1]
heights <- rep(magn,nrows)
order <- matrix( c(1,2,3), nrows, ncols, byrow=TRUE )

pdf( "fig/gpp_loss.pdf", width=sum(widths), height = sum(heights) )
  
  par(las=1)

  panel <- layout(
          order,
          widths=widths,
          heights=heights,
          TRUE
          )
  # layout.show( panel )

  ##------------------------------------------------------------------------
  ## map
  ##------------------------------------------------------------------------
  out <- gpp_s0 - gpp_s1b
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
  
  lev = c( 0, 400, 10 )
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=max(out, na.rm=TRUE), minval=min(out, na.rm=TRUE) )
  
  par( mar=c(4,3,3,1),xaxs="i", yaxs="i",las=1)
  image(
          lon, lat,
          out,
          ylim=c(-60,85),
          yaxt="n", xaxt="n",
          col=out.mycolorbar$colors, breaks=out.mycolorbar$margins,
          xlab="", ylab=""
          )
  map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
  
  axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
  axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 4, at=lat.labels, lab=F, lwd=1.5 )
  axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
  axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )
  
  axis( 3, at=lon.labels, lab=F, lwd=1.5 )
  axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )
  
  mtext( "a)", font=2, adj = 0, line = 0.5, cex = 1.2 )

  ## Color key
  par( mar=c(4,3,3,1),xaxs="i", yaxs="i",las=1)
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1900 )

  ##------------------------------------------------------------------------
  ## time series
  ##------------------------------------------------------------------------
  par( mar=c(4,4.1,3,1))
  with( df, plot( year, gpp_s1b, type="l", ylim=c(100,160), lty=2, xlab="Year", ylab=expression( paste("Global GPP (Pg C yr"^-1, ")" ) ) ) )
  with( df,  polygon( c( year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), border = NA, col=rgb(0,0,0,0.3) ) )

  with( df, lines( year, gpp_s0 ) )
  with( df, lines( year, gpp_modis, col="red" ) )
  with( df, lines( year, gpp_vpm, col="magenta" ) )
  with( df, lines( year, gpp_bess, col="blue" ) )
  with( df, lines( year, gpp_mte, col="green" ) )

  legend("topleft", c( "P-model", "P-model, corrected (IV)", "MODIS", "VPM", "BESS", "MTE" ), bty = "n", lty = c(1,2,1,1,1,1), col=c("black", "black", "red", "magenta", "blue", "green") )
  mtext( "b)", font=2, adj = 0, line = 0.5, cex = 1.2 )

dev.off()