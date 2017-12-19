library(dplyr)
library(ncdf4)

nc <- nc_open( paste0( myhome, "data/trendy/v5/LPX-Bern/S2/LPX_S2_gpp.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$longitude$vals
lat <- nc$dim$latitude$vals
nc_close(nc)

## get time from a different file
nc <- nc_open( paste0( myhome, "data/trendy/v5/CLM/S2/CLM4.5_S2_gpp.nc" ) )
time_units <- nc$dim$time$units
time <- nc$dim$time$vals
nc_close(nc)

gpp <- array( gpp, dim = c(360,180,1,1872))

cdf.write( gpp, "gpp",
           lon, lat,
           filnam = paste0( myhome, "data/trendy/v5/LPX-Bern/S2/LPX_S2_gpp_NICE.nc" ),
           nvars = 1,
           time = time,
           make.zdim = TRUE,
           z_dim = 1,
           make.tdim = TRUE,
           units_time = time_units,
           units_var1 = "kg C m-2 s-1",
           glob_hist = "created by soilm_global/preproc_lpx.R"
)
