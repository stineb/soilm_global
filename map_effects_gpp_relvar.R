library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

## file name for figure
filn <- "fig/map_gpp_relvar_diff.pdf"

##------------------------------------------------------------------------
## GPP interannual relative variance change
##------------------------------------------------------------------------
fil_s0  <- "gpp_pmodel_s0_RELVAR.nc"
fil_s1  <- "gpp_pmodel_s1_RELVAR.nc"
fil_s1a <- "gpp_pmodel_s1a_RELVAR.nc"
fil_s1b <- "gpp_pmodel_s1b_RELVAR.nc"
fil_s1c <- "gpp_pmodel_s1c_RELVAR.nc"

dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

## S0
nc <- nc_open( paste0( dir, fil_s0 ) )
gpp_s0 <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

## S1
nc <- nc_open( paste0( dir, fil_s1 ) )
gpp_s1 <- ncvar_get( nc, varid="gpp" )
nc_close(nc)

## S1a
nc <- nc_open( paste0( dir, fil_s1a ) )
gpp_s1a <- ncvar_get( nc, varid="gpp" )
nc_close(nc)

## S1b
nc <- nc_open( paste0( dir, fil_s1b ) )
gpp_s1b <- ncvar_get( nc, varid="gpp" )
nc_close(nc)

## S1c
nc <- nc_open( paste0( dir, fil_s1c ) )
gpp_s1c <- ncvar_get( nc, varid="gpp" )
nc_close(nc)


##-----------------------------------------------------
## Plot relative variance in S1
##-----------------------------------------------------
plot_map( gpp_s1, lev=c( 0, 40, 10 ),
      toplefttext=expression(paste("GPP relative variance")),
      toprighttext=expression(paste("fraction")),
      maxval = 200, color=cols
      )

##-----------------------------------------------------
## Change in relative variance
##-----------------------------------------------------
  arr = gpp_s1 / gpp_s0
  toplefttext = expression(paste("GPP amplification of relative variance"))
  toprighttext = expression(paste("fraction"))
  minval = NA
  maxval = 35
  color = c( "royalblue4", "wheat", "tomato2", "tomato4" )
  lev=c( 0, 4, 10 )

  ## half degree resolution
  lon <- seq(-179.75, 179.75, 0.5)
  lat <- seq(-89.75, 89.75, 0.5)

  magn <- 4
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[1] <- 0.15*widths[1]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  if (!is.na(filn)) pdf( filn, width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

    ## Color key
    color <- c( "royalblue4", "wheat", "tomato2", "tomato4" )
    lev <- c( 0, 4, 10 )
    maxval = 35
    minval = NA
    par( mar=c(2.5,3,1,1),xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=maxval, minval=minval )

    par( mar=c(2.5,2.5,1,1),xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
    image(
            lon, lat, 
            arr,
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

    mtext( expression(paste("GPP amplification of relative variance")), line=1, adj=0 )
    mtext( expression(paste("fraction")), line=1, adj=1 )

    rect( -179, -50, -100, 7.5, border = NA, col="white" )
    
    ##-----------------------------------------------------
    ## Inset: Empirical cumulative distribution function of the amplification factor
    ##-----------------------------------------------------
    vec <- c( gpp_s1 / gpp_s0 )
    vec <- vec[!is.na(vec)]
    ecdf_ampl <- ecdf( vec )

    vec <- c( gpp_s1a / gpp_s0 )
    vec <- vec[!is.na(vec)]
    ecdf_ampl_a <- ecdf( vec )

    vec <- c( gpp_s1b / gpp_s0 )
    vec <- vec[!is.na(vec)]
    ecdf_ampl_b <- ecdf( vec )

    vec <- c( gpp_s1c / gpp_s0 )
    vec <- vec[!is.na(vec)]
    ecdf_ampl_c <- ecdf( vec )

    ## Inset 1
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v_orig <- v
    v <- c( v[1]+0.07, v[1]+0.2*v[2], v[3]+0.10*v[4], v[3]+0.32*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.2,0) )
    xlim <- c(0.75,20)
    ylim <- c(0.001, 1)
    plot( xlim, ylim, type="n", xlim=xlim, log="xy", ylim=ylim, xlab = "", ylab = "", bg="white", cex.axis=0.7, tck=-0.03 )
    mtext( "amplification", side=1, line=1, adj=0.5, cex = 0.7 )
    mtext( "ECDF", side=2, line=1.7, adj=0.5, cex = 0.7, las=0 )
    rect( xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col="white" )
    curve( 1.0 - ecdf_ampl(x), from=xlim[1], to=xlim[2], col="red", add=TRUE  )
    polygon( c( seq(xlim[1], xlim[2], by=0.1 ), rev( seq(xlim[1], xlim[2], by=0.1 ) ) ),  c( 1 - ecdf_ampl_a( seq(xlim[1], xlim[2], by=0.1 ) ), rev( 1- ecdf_ampl_c( seq(xlim[1], xlim[2], by=0.1 ) ) ) ), border = NA, col = rgb(1,0,0,0.4) )
    abline( v=1, lty=3 )
    box()

  if (!is.na(filn)) dev.off()

## save vector amplification factor
save( vec, file="data/ampl_relvar.Rdata" )
