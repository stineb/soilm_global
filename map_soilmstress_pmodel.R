library(ncdf4)
library(fields)
library(sp)
library(maptools)
library(dplyr)

##------------------------------------------------------------------------
## GPP loss
##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_mean_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_mean_s1_fapar3g_global.nc"

  dir <- "/Users/benjaminstocker/data/pmodel_fortran_output/"


  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  nc <- nc_open( paste0( dir, fil_s1 ) )
  gpp_s1 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  ##-----------------------------------------------------
  ## absolute GPP loss (gC m-2)
  ##-----------------------------------------------------
    magn <- 4
    ncols <- 2
    nrows <- 1
    widths <- rep(1.6*magn,ncols)
    widths[2] <- 0.18*widths[1]
    heights <- rep(magn,nrows)
    order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)

    a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
    b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

    # pdf( "fig/map_gpploss_abs_pmodel.pdf", width=sum(widths), height=sum(heights) )

      panel <- layout(
                order,
                widths=widths,
                heights=heights,
                TRUE
                )
      # layout.show( panel )

      ## Color key
      # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      color <- c( "wheat", "tomato2", "tomato4" )
      lev <- c( 0, 600, 10 )
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=1100 )

      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      image(
              seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
              gpp_s0 - gpp_s1,
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

      mtext( expression(paste("GPP loss")), line=1, adj=0, font=2 )
      mtext( expression(paste("gC m"^{-2}, "yr"^{-1})), line=1, adj=1 )

      ## Color key
      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

    # dev.off()


  ##-----------------------------------------------------
  ## relative GPP loss (percent)
  ##-----------------------------------------------------
    magn <- 4
    ncols <- 2
    nrows <- 1
    widths <- rep(1.6*magn,ncols)
    widths[2] <- 0.18*widths[1]
    heights <- rep(magn,nrows)
    order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)

    a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
    b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

    # pdf( "fig/map_gpploss_abs_pmodel.pdf", width=sum(widths), height=sum(heights) )

      panel <- layout(
                order,
                widths=widths,
                heights=heights,
                TRUE
                )
      # layout.show( panel )

      ## Color key
      # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      color <- c( "wheat", "tomato2", "tomato4" )
      lev <- c( 0, 50, 10 )
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=100 )

      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      image(
              seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
              (1-(gpp_s1 / gpp_s0))*100,
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

      mtext( expression(paste("GPP loss")), line=1, adj=0, font=2 )
      mtext( expression(paste("%")), line=1, adj=1 )

      ## Color key
      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

    # dev.off()


##------------------------------------------------------------------------
## GPP interannual relative variance change
##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_relvar_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_relvar_s1_fapar3g_global.nc"

  dir <- "/Users/benjaminstocker/data/pmodel_fortran_output/"


  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  nc <- nc_open( paste0( dir, fil_s1 ) )
  gpp_s1 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  ##-----------------------------------------------------
  ## Change in relative variance
  ##-----------------------------------------------------
    magn <- 4
    ncols <- 2
    nrows <- 1
    widths <- rep(1.6*magn,ncols)
    widths[2] <- 0.18*widths[1]
    heights <- rep(magn,nrows)
    order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)

    a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
    b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

    # pdf( "fig/map_gpploss_abs_pmodel.pdf", width=sum(widths), height=sum(heights) )

      panel <- layout(
                order,
                widths=widths,
                heights=heights,
                TRUE
                )
      # layout.show( panel )

      ## Color key
      # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      color <- c( "royalblue4", "wheat", "tomato2", "tomato4" )
      lev <- c( 0, 4, 10 )
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=35 )

      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      image(
              seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
              gpp_s1 / gpp_s0,
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

      mtext( expression(paste("GPP change in relative variance")), line=1, adj=0, font=2 )
      mtext( expression(paste("fraction")), line=1, adj=1 )

      ## Color key
      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )


##------------------------------------------------------------------------
## GPP interannual (absolute) variance change
##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_var_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_var_s1_fapar3g_global.nc"

  dir <- "/Users/benjaminstocker/data/pmodel_fortran_output/"


  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  nc <- nc_open( paste0( dir, fil_s1 ) )
  gpp_s1 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  ##-----------------------------------------------------
  ## Difference in absolute variance
  ##-----------------------------------------------------
    magn <- 4
    ncols <- 2
    nrows <- 1
    widths <- rep(1.6*magn,ncols)
    widths[2] <- 0.18*widths[1]
    heights <- rep(magn,nrows)
    order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)

    a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
    b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

    # pdf( "fig/map_gpploss_abs_pmodel.pdf", width=sum(widths), height=sum(heights) )

      panel <- layout(
                order,
                widths=widths,
                heights=heights,
                TRUE
                )
      # layout.show( panel )

      ## Color key
      # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      color <- c( "royalblue4", "royalblue2", "wheat", "tomato2", "tomato4" )
      lev <- c( -10, 10, 10 )
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=50, minval=-50 )

      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      image(
              seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
              (gpp_s1 - gpp_s0)*1e-3,
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

      mtext( expression(paste("GPP difference in absolute variance")), line=1, adj=0, font=2 )
      mtext( expression(paste("gC m"^{-2}, "yr"^{-1})), line=1, adj=1 )

      ## Color key
      par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
      out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

