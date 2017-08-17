library(ncdf4)
library(fields)
library(sp)
library(maptools)
library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

## get GPP
if (!exists("gpp")){
  print("getting gpp ...")
  nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_GPP_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  gpp <- ncvar_get( nc, varid="GPP" )
  nc_close( nc )
}

## get annual mean alpha (AET/PET)
if (!exists("aalpha")){
  print("getting gpp ...")
  nc <- nc_open( "data/alpha_Pmod_SPLASH.nc" )
  aalpha <- ncvar_get( nc, varid="ALPHA" )
  nc_close( nc )
}



## get %GPP loss as a function of AET/PET from Fig. 9a in Stocker et al. in prep.
calc_gpploss <- function( alpha ){
  out <- ( 52.3 - 54.1 * alpha )
  return( out )
}
gpploss_frac <- apply( aalpha, c(1,2), FUN = calc_gpploss )

gpploss_frac[ which( gpploss_frac>100 ) ] <- 100
gpploss_frac[ which( gpploss_frac<0 ) ]   <- 0

## get mean annual GPP
print("get mean annual GPP...")
agpp <- apply( gpp, c(1,2), FUN = function(x) sum( x, na.rm=TRUE )/12 )

## get absolute GPP loss
print("get annual GPP loss ...")
gpploss <- agpp * gpploss_frac * 1e-3


##------------------------
## %GPP loss, map
##------------------------
  magn <- 4
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[2] <- 0.25*widths[1]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  pdf( "fig/map_gpploss_frac.pdf", width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

    ## Color key
    # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    color <- c( "wheat", "tomato" )
    lev <- c( 0, 50, 10 )
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=100 )

    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    image(
            seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
            gpploss_frac,
            ylim=c(-60,85), 
            # zlim=range(lev), 
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

    ## Color key
    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

  dev.off()


##------------------------
## absolute GPP loss, map
##------------------------
  magn <- 4
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[2] <- 0.25*widths[1]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  pdf( "fig/map_gpploss_abs.pdf", width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

    ## Color key
    # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    color <- c( "wheat", "tomato" )
    lev <- c( 0, 100, 10 )
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=1000 )

    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    image(
            seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
            gpploss,
            ylim=c(-60,85), 
            # zlim=range(lev), 
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

    ## Color key
    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

  dev.off()  
