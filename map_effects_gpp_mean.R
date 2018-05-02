library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP loss
##------------------------------------------------------------------------
  fil_s0  <- "gpp_pmodel_s0_MEAN.nc"
  fil_s1  <- "gpp_pmodel_s1_MEAN.nc"
  fil_s1b <- "gpp_pmodel_s1b_MEAN.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  # ## S1
  # nc <- nc_open( paste0( dir, fil_s1 ) )
  # gpp_s1 <- ncvar_get( nc, varid="gpp" )
  # nc_close(nc)

  ## S1b
  nc <- nc_open( paste0( dir, fil_s1b ) )
  gpp_s1b <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)


  ##-----------------------------------------------------
  ## absolute GPP (gC m-2 yr-1)
  ##-----------------------------------------------------
  cols <- colorRampPalette( brewer.pal(9,"YlOrRd"))(10)

  plot_map( gpp_s0, lev = c( 0, 4500, 10 ), 
            toplefttext=expression(paste("GPP, no soil moisture limitation")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=4600, color=cols
            , file = "fig/map_gpp_nolimit.pdf"
            )

  plot_map( gpp_s1b, lev = c( 0, 4500, 10 ), 
            toplefttext=expression(paste("GPP, with soil moisture limitation")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=4600, color=cols
            , file = "fig/map_gpp_limit.pdf"
            )

  ##-----------------------------------------------------
  ## absolute GPP loss (gC m-2 yr-1)
  ##-----------------------------------------------------
  plot_map( gpp_s0 - gpp_s1b, lev = c( 0, 600, 10 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=1900, color=cols
            , file = "fig/map_gpp_loss_abs.pdf"
            )

  ##-----------------------------------------------------
  ## relative GPP loss (percent)
  ##-----------------------------------------------------
  plot_map( (1-(gpp_s1b / gpp_s0))*100, lev=c( 0, 70, 14 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("%")),
            color=cols, maxval = 100
            , file = "fig/map_gpp_loss_rel.pdf"
            )


##------------------------------------------------------------------------
## Plot for paper
##------------------------------------------------------------------------
require( ncdf4, quietly = TRUE )
require( fields, quietly = TRUE )
require( sp, quietly = TRUE )
require( maptools, quietly = TRUE )
require( dplyr, quietly = TRUE )  

## load global GPP time series, prepared by plot_effects_gpp_tseries.R
load("data/gpp_glob_tseries.Rdata")

magn <- 4
ncols <- 3
nrows <- 1
widths <- rep(1.3*magn,ncols)
widths[2] <- 0.15*widths[1]
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

  ## map
  out <- gpp_s0 - gpp_s1b
  lon <- seq(-179.75, 179.75, 0.5)
  lat <- seq(-89.75, 89.75, 0.5)
  
  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)
  
  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )
  
  color <- c( "wheat", "tomato2", "tomato4" )
  
  lev = c( 0, 500, 10 )
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=max(out, na.rm=TRUE), minval=min(out, na.rm=TRUE) )
  
  par( mar=c(3,3,3,1),xaxs="i", yaxs="i",las=1)
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
  
  mtext( "a)", font=2, adj = 0, line = 0.5, cex = 1 )

  ## Color key
  par( mar=c(3,3,3,1),xaxs="i", yaxs="i",las=1)
  out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1900 )
  
  ## time series
  par( mar=c(3,4,3,1))
  with( df, plot( year, gpp_s1b, type="l", ylim=c(100,160), lty=2, xlab="year", ylab=expression( paste("Global GPP (PgC yr"^-1, ")" ) ) ) )
  with( df,  polygon( c( year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), border = NA, col=rgb(0,0,0,0.3) ) )

  with( df, lines( year, gpp_s0 ) )
  with( df, lines( year, gpp_modis, col="red" ) )
  with( df, lines( year, gpp_vpm, col="magenta" ) )
  with( df, lines( year, gpp_bess, col="blue" ) )
  with( df, lines( year, gpp_mte, col="green" ) )
  axis(4, labels = FALSE)

  legend("topleft", c( "P-model, s0", "P-model, s1b", "MODIS", "VPM", "BESS", "MTE" ), bty = "n", lty = c(1,2,1,1,1,1), col=c("black", "black", "red", "magenta", "blue", "green") )
  mtext( "b)", font=2, adj = 0, line = 0.5, cex = 1 )

dev.off()

