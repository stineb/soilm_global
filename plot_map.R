plot_map <- function( arr, lev, file=NA, positive=TRUE, toplefttext=NA, toprighttext=NA, minval=NA, maxval=NA, color=NA ){

  require( ncdf4, quietly = TRUE )
  require( fields, quietly = TRUE )
  require( sp, quietly = TRUE )
  require( maptools, quietly = TRUE )
  require( dplyr, quietly = TRUE )  

  if ( dim(arr)[1]==720 && dim(arr)[2]==360 ){

    ## half degree resolution
    lon <- seq(-179.75, 179.75, 0.5)
    lat <- seq(-89.75, 89.75, 0.5)

  } else if ( dim(arr)[1]==360 && dim(arr)[2]==180 ){

    ## one degree resolution
    lon <- seq(-179.5, 179.5, 1.0 )
    lat <- seq(-89.5, 89.5, 1.0 ) 

  }    

  magn <- 4
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[2] <- 0.15*widths[1]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  if (!is.na(file)) pdf( file, width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

    ## Color key
    if (is.na(color)){
      if (positive){
        color <- c( "wheat", "tomato2", "tomato4" )
      } else {
        color <- c( "royalblue4", "royalblue2", "wheat", "tomato2", "tomato4" )      
      }
    }

    # lev <- c( 0, 0.5, 0.7, 0.8, 0.9, 0.95, 1, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2, 999 )
    # lev <- seq(-2,2,0.2)
    # lev <- c( 0, 0.2, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.5, 2, 3, 999 )

    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=maxval, minval=minval )

    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
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

    if (!is.na(toplefttext)) mtext( toplefttext, line=1, adj=0 )
    if (!is.na(toprighttext)) mtext( toprighttext, line=1, adj=1 )

    ## Color key
    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

  if (!is.na(file)) dev.off()  

}