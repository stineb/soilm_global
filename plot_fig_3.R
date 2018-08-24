library( RColorBrewer )
library( fields, quietly = TRUE )
library( sp, quietly = TRUE )
library( maptools, quietly = TRUE )
library( dplyr, quietly = TRUE )  

source("../utilities/add_alpha.R")

load("data/relvar.Rdata")
load("data/ampl_agg_jung.Rdata")

## file name for figure
filn <- "fig/fig_3.pdf"

##-----------------------------------------------------
## Panel layout
##-----------------------------------------------------
  magn <- 4
  ncols <- 3
  nrows <- 1
  widths <- rep(1.4*magn,ncols)
  widths[2] <- 0.17*widths[1]
  widths[1] <- 0.80*widths[3]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2,3), nrows, ncols, byrow=TRUE )

  parinit <- par( no.readonly=TRUE )
  if (!is.na(filn)) pdf( filn, width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

  ##-----------------------------------------------------
  ## Jung plot
  ##-----------------------------------------------------
    vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )
    
    par( parinit, mar=c(4,4,3,1), xaxs="i", yaxs="i",las=1, mgp=c(3,1,0) )
    # par(las=1, mar=c(4,4,1,1), new=FALSE, fig=c(0, 1, 0, 1) )
    # with( ampl_agg, plot( resnr, relvar_median, type="l", col="black", lwd=2, ylim=c(0,10), xaxt = "n", xlab="Spatial resolution (degrees)", ylab="Amplification" ) )
    with( ampl_agg, plot( resnr, relvar_median, type="l", col="black", lwd=2, ylim=c(0,10), xaxt = "n", xlab=bquote( "Spatial resolution "(degree) ), ylab="Amplification" ) )
    axis( 1, at=seq(length(vec_res)), labels=as.character( vec_res ) )
    with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q01, relvar_q99), col=add_alpha("tomato2", 0.25), border = NA ) )
    with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q05, relvar_q95), col=add_alpha("tomato2", 0.25), border = NA ) )
    with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q10, relvar_q90), col=add_alpha("tomato2", 0.25), border = NA ) )
    with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q25, relvar_q75), col=add_alpha("tomato2", 0.25), border = NA ) )
    abline( h = 1.0, lty=3 )
    mtext( "a)", font=2, adj = 0, line = 0.5, cex = 1 )

    
  ##-----------------------------------------------------
  ## Map change in relative variance
  ##-----------------------------------------------------
    arr = gpp_s1b / gpp_s0

    # toplefttext = expression(paste("GPP amplification of relative variance"))
    # toprighttext = expression(paste("fraction"))
    minval = NA
    maxval = 35
    color = c( "royalblue4", "wheat", "tomato2", "tomato4" )
    lev=c( 0, 4, 10 )

    ## half degree resolution
    lon <- seq(-179.75, 179.75, 0.5)
    lat <- seq(-89.75, 89.75, 0.5)

    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)

    a <- sapply( lat.labels, function(x) if (x>0) {bquote(.(x)*degree ~N)} else if (x==0) {bquote(.(x)*degree)} else {bquote(.(-x)*degree ~S)} )
    b <- sapply( lon.labels, function(x) if (x>0) {bquote(.(x)*degree ~E)} else if (x==0) {bquote(.(x)*degree)} else {bquote(.(-x)*degree ~W)})

    color <- c( "royalblue4", "wheat", "tomato2", "tomato4" )
    lev <- c( 0, 4, 10 )
    maxval = 35
    minval = NA
    par( mar=c(4,3,3,1),xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=maxval, minval=minval )

    par( mar=c(4,2.5,3,1), xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
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

    # mtext( expression(paste("GPP amplification of relative variance")), line=1, adj=0 )
    # mtext( expression(paste("fraction")), line=1, adj=1 )

    rect( -179, -50, -100, 7.5, border = NA, col="white" )

    mtext( "b)", font=2, adj = 0, line = 0.5, cex = 1 )

  ##-----------------------------------------------------
  ## Inset: Empirical cumulative distribution function of the amplification factor
  ##-----------------------------------------------------
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
    v <- c( v[1]+0.035, v[1]+0.12*v[2], v[3]+0.10*v[4], v[3]+0.36*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.2,0) )
    xlim <- c(0.75,25)
    ylim <- c(0.001, 1)
    plot( xlim, ylim, type="n", xlim=xlim, log="xy", ylim=ylim, xlab = "", ylab = "", bg="white", cex.axis=0.7, tck=-0.03 )
    mtext( "Amplification", side=1, line=1, adj=0.5, cex = 0.7 )
    mtext( "ECDF", side=2, line=1.7, adj=0.5, cex = 0.7, las=0 )
    rect( xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col="white" )
    curve( 1.0 - ecdf_ampl_b(x), from=xlim[1], to=xlim[2], col="red", add=TRUE  )
    polygon( c( seq(xlim[1], xlim[2], by=0.1 ), rev( seq(xlim[1], xlim[2], by=0.1 ) ) ),  c( 1 - ecdf_ampl_a( seq(xlim[1], xlim[2], by=0.1 ) ), rev( 1- ecdf_ampl_c( seq(xlim[1], xlim[2], by=0.1 ) ) ) ), border = NA, col = rgb(1,0,0,0.4) )
    abline( v=1, lty=3 )
    box()


  if (!is.na(filn)) dev.off()

