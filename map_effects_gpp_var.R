library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP interannual (absolute) variance change
##------------------------------------------------------------------------
fil_s0 <- "gpp_pmodel_s0_VAR.nc"
fil_s1 <- "gpp_pmodel_s1_VAR.nc"

dir <- paste0( myhome, "/data/pmodel_fortran_output/")

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

##-----------------------------------------------------
## Plot absolute variance in S1
##-----------------------------------------------------
plot_map( gpp_s1*1e-3, lev=c( 0, 40, 10 ), 
      toplefttext=expression(paste("GPP aboslute variance")), 
      toprighttext=expression(paste("gC m"^{-2}, "yr"^{-1})),
      maxval = 200, color=cols
      , file = "fig/map_gpp_var_s1.pdf" 
      )

##-----------------------------------------------------
## Difference in absolute variance
##-----------------------------------------------------
plot_map( (gpp_s1 - gpp_s0)*1e-3, lev=c( -10, 10, 10 ), 
          toplefttext=expression(paste("GPP change in absolute variance")), 
          toprighttext=expression(paste("kgC m"^{-2}, "yr"^{-1})),
          positive = FALSE, maxval=60, minval=-60
          , file = "fig/map_gpp_var_diff.pdf"
          )



