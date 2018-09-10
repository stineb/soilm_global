library(ncdf4)
library(dplyr)
library(RColorBrewer)
source("plot_map.R")

get_trend <- function( vec ){
  if (is.na(vec[1])){
    slope <- NA
    signif <- FALSE
  } else {
    df <- tibble( x = seq(length(vec)), y = vec )
    linmod <- lm( y ~ x, data = df )
    slope <- coef(linmod)[2]
    ci <- confint(linmod)[2,]
    # if (ci[1] < 0.0 && 0.0 < ci[2]) slope <- NA
    if (ci[1] < 0.0 && 0.0 < ci[2]) {
      signif <- FALSE
    } else {
      signif <- TRUE
    }
  }
  return( c(slope, signif) )
}

##------------------------------------------------------------------------
## TREND in relative GPP reduction, across globe
##------------------------------------------------------------------------
  fil_s0  <- "gpp_pmodel_s0_ANN.nc"
  fil_s1b <- "gpp_pmodel_s1b_ANN.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)

  ## S1b
  nc <- nc_open( paste0( dir, fil_s1b ) )
  gpp_s1b <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)

  ## get relative reduction in annual GPP due to soil moisture
  gpp_reldiff <- -100 * ( 1.0 - gpp_s1b/gpp_s0 )
  gpp_absdiff <- gpp_s1b - gpp_s0
  
  ## get linear trend in each gridcell
  # trend_abs <- apply( gpp_absdiff, c(1,2), FUN = get_trend )
  trend_rel <- apply( gpp_reldiff, c(1,2), FUN = get_trend )

  # ## plot absolute difference
  # plot_map( trend_abs[1,,], lev = c(-5, 5, 10), positive = FALSE, minval = -20, maxval = 40 )
  
  ## plot relative difference (positive value meaning that the relative difference due to soil moisture is getting less negative = less of a reduction in GPP in relative terms)
  plot_map( trend_rel[1,,], lev = c(-0.5, 0.5, 10), positive = FALSE, minval = -1.4, maxval = 2.0 ) # , stippling =  trend_rel[2,,]
  
##------------------------------------------------------------------------
## Get mean annual GPP fields
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

  ## S1
  nc <- nc_open( paste0( dir, fil_s1 ) )
  gpp_s1 <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)

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
  plot_map( gpp_s0 - gpp_s1b, lev = c( 0, 400, 10 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=1900
            , file = "fig/map_gpp_loss_abs.pdf"
            )

  ##-----------------------------------------------------
  ## relative GPP loss (percent)
  ##-----------------------------------------------------
  plot_map( (1-(gpp_s1b / gpp_s0))*100, lev=c( 0, 70, 14 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("%")),
            maxval = 100
            , file = "fig/map_gpp_loss_rel.pdf"
            )

  save( gpp_s0, gpp_s1b, trend_rel, file="data/gpp_loss.Rdata" )

