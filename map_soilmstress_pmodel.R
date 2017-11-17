library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP loss
##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_mean_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_mean_s1_fapar3g_global.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/")

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
  plot_map( (1-(gpp_s1 / gpp_s0))*100, lev=c( 0, 50, 10 ), 
            toplefttext=expression(paste("GPP loss")), 
            toprighttext=expression(paste("%")),
            color=cols, maxval = 100
            )


##------------------------------------------------------------------------
## GPP interannual relative variance change
##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_relvar_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_relvar_s1_fapar3g_global.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/")

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
  plot_map( gpp_s1 / gpp_s0, lev=c( 0, 4, 10 ), 
        toplefttext=expression(paste("GPP change in relative variance")), 
        toprighttext=expression(paste("fraction")),
        color = c( "royalblue4", "wheat", "tomato2", "tomato4" ),
        maxval = 35
        )


  ##------------------------------------------------------------------------
  ## GPP interannual (absolute) variance change
  ##------------------------------------------------------------------------
  fil_s0 <- "pmodel_gpp_var_s0_fapar3g_global.nc"
  fil_s1 <- "pmodel_gpp_var_s1_fapar3g_global.nc"
  
  dir <- paste0( myhome, "/data/pmodel_fortran_output/")
  
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
  ## Plot absolute variance in S1
  ##-----------------------------------------------------
  plot_map( gpp_s1*1e-3, lev=c( 0, 40, 10 ), 
        toplefttext=expression(paste("GPP aboslute variance")), 
        toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
        maxval = 200, color=cols
        )
  
  ##-----------------------------------------------------
  ## Difference in absolute variance
  ##-----------------------------------------------------
  plot_map( (gpp_s1 - gpp_s0)*1e-3, lev=c( -10, 10, 10 ), 
            toplefttext=expression(paste("GPP difference in absolute variance")), 
            toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
            positive = FALSE, maxval=60, minval=-60
            )

