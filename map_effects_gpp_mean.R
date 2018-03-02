library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP loss
##------------------------------------------------------------------------
  fil_s0 <- "gpp_pmodel_s0_MEAN.nc"
  fil_s1 <- "gpp_pmodel_s1_MEAN.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/")

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


  ##-----------------------------------------------------
  ## absolute GPP (gC m-2 yr-1)
  ##-----------------------------------------------------
  cols <- colorRampPalette( brewer.pal(9,"YlOrRd"))(10)

  plot_map( gpp_s0, lev = c( 0, 4500, 10 ), 
            toplefttext=expression(paste("GPP, no soil moisture limitation")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=4600, color=cols
            )

  plot_map( gpp_s1, lev = c( 0, 4500, 10 ), 
            toplefttext=expression(paste("GPP, with soil moisture limitation")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=4600, color=cols
            )

  ##-----------------------------------------------------
  ## absolute GPP loss (gC m-2 yr-1)
  ##-----------------------------------------------------
  plot_map( gpp_s0 - gpp_s1, lev = c( 0, 600, 10 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            maxval=1900, color=cols
            )

  ##-----------------------------------------------------
  ## relative GPP loss (percent)
  ##-----------------------------------------------------
  plot_map( (1-(gpp_s1 / gpp_s0))*100, lev=c( 0, 70, 14 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("%")),
            color=cols, maxval = 100
            )
