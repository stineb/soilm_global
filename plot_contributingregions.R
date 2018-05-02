library(ncdf4)
library(abind)

source("../utilities/plot_map.R")


get_stocker_f <- function( eff, anom, isabs=FALSE ){
  ##---------------------------------------------------
  ## requires as input a 3D effay with lon x lat x time
  ## and values being annual detrended anomalies 
  ## index quantifies the degree to which each gridcell contributes to the global signal
  ##---------------------------------------------------

  ## use gridcell total not per unit area
  if (isabs==FALSE){
    source( "integrate_gridcell.R" )
    
    eff_abs <- integrate_gridcell( eff, global=FALSE, overwrite=TRUE )
    anom_abs <- integrate_gridcell( anom, global=FALSE, overwrite=TRUE )

    eff_glob <- apply( eff_abs, c(3), FUN=sum, na.rm=TRUE )
    anom_glob <- apply( anom_abs, c(3), FUN=sum, na.rm=TRUE )
  } else {
    
    eff_abs <- eff
    eff_glob <- apply( eff_abs, c(3), FUN=sum, na.rm=TRUE )    
    anom_glob <- apply( anom_abs, c(3), FUN=sum, na.rm=TRUE )    
  }

  stocker_f <- eff[,,1]
  stocker_f[] <- NA
  for (ilon in seq(dim(eff)[1])){
    for (ilat in seq(dim(eff)[2])){
      if (!is.na(eff[ilon,ilat,1])){
        stocker_f[ilon,ilat] <- sum( eff_abs[ilon,ilat,] * abs( anom_glob ) / anom_glob ) / sum( abs( anom_glob ) )
      }
    }
  }
  return( stocker_f )
}


filpath_detr <- c(  paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s0_DETR.nc"), 
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_DETR.nc")
                    )

modl <- c( "Pmodel_S0", "Pmodel_S1a", "Pmodel_S1b", "Pmodel_S1c")

detr <- list()
for (idx in seq(length(modl))){

  ## Read detrended GPP for IAV
  if (file.exists(filpath_detr[idx])){

    ## read file
    nc <- nc_open( filpath_detr[idx] )
    detr[[ modl[idx] ]] <- try( ncvar_get( nc, varid="gpp" ) )
    nc_close(nc)

  }
  
}

## get difference
detr[[ "diffa" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1a" ]]
detr[[ "diffb" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1b" ]]
detr[[ "diffc" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1c" ]]

stocker_fb <- get_stocker_f( detr$diffb, detr$Pmodel_S1b, isabs=FALSE )

# for (it in 1:5){
#   ahlstroem_fb <- get_stocker_f( detr$Pmodel_S1b[,,((it-1)*7+1):((it-1)*7+7)], detr$Pmodel_S1b[,,((it-1)*7+1):((it-1)*7+7)], isabs=FALSE )
#   plot_map( ahlstroem_fb*1e4, lev=seq(-1,1,0.2), positive=FALSE, maxval=30, minval=-30, file=paste0("fig/map_ahlstroem", as.character(it), ".pdf") ) #  
# }
# ahlstroem_fb <- get_stocker_f( detr$Pmodel_S1b, detr$Pmodel_S1b, isabs=FALSE )
# plot_map( ahlstroem_fb*1e4, lev=seq(-1,1,0.2), positive=FALSE, maxval=30, minval=-30, file=paste0("fig/map_ahlstroem.pdf") ) #  


# ahlstroem_fb <- get_stocker_f( detr$Pmodel_S0, detr$Pmodel_S0, isabs=FALSE )
# plot_map( ahlstroem_fb*1e4, lev=seq(-1,1,0.2), positive=FALSE, maxval=30, minval=-30, file=paste0("fig/map_ahlstroem_s0.pdf") ) #  

##-----------------------------------------------------
## Plot without inset
##-----------------------------------------------------
plot_map( stocker_fb*1e4, lev=seq(-0.5,0.5,0.1), positive=FALSE, maxval=30, minval=-30 )

##-----------------------------------------------------
## Plot with inset
##-----------------------------------------------------
filn <- "fig/map_stocker_gpploss.pdf"
  arr = stocker_fb*1e4
  toplefttext = expression(paste(""))
  toprighttext = expression(paste(""))
  minval = -30
  maxval = 30
  color <- c( "royalblue4", "royalblue2", "wheat", "tomato2", "tomato4" )  
  lev = seq(-0.5,0.5,0.1)

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
    ## Inset 1
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v_orig <- v
    v <- c( v[1]+0.03, v[1]+0.2*v[2], v[3]+0.10*v[4], v[3]+0.32*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.5,0) )
    hist( arr, breaks=50, main="", freq = FALSE, xlim=c(-0.75,0.75), cex.axis=0.7, axes=FALSE, col="grey70" )
    axis( 1, cex.axis=0.7 )

    # plot( xlim, ylim, type="n", xlim=xlim, log="xy", ylim=ylim, xlab = "", ylab = "", bg="white", cex.axis=0.7, tck=-0.03 )
    # mtext( "amplification", side=1, line=1, adj=0.5, cex = 0.7 )
    # mtext( "ECDF", side=2, line=1.7, adj=0.5, cex = 0.7, las=0 )
    # rect( xlim[1], ylim[1], xlim[2], ylim[2], border = NA, col="white" )
    # curve( 1.0 - ecdf_ampl(x), from=xlim[1], to=xlim[2], col="red", add=TRUE  )
    # polygon( c( seq(xlim[1], xlim[2], by=0.1 ), rev( seq(xlim[1], xlim[2], by=0.1 ) ) ),  c( 1 - ecdf_ampl_a( seq(xlim[1], xlim[2], by=0.1 ) ), rev( 1- ecdf_ampl_c( seq(xlim[1], xlim[2], by=0.1 ) ) ) ), border = NA, col = rgb(1,0,0,0.4) )
    # abline( v=1, lty=3 )
    # box()

  if (!is.na(filn)) dev.off()

