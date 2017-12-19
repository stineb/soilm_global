library(dplyr)
library(ncdf4)

source( paste0( myhome, "/utilities/get_days_since.R" ) )

nc <- nc_open( paste0( myhome, "data/trendy/v5/SDGVM/S2/SDGVM_S2_gpp.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
nc_close(nc)

gpp <- array( gpp, dim = c(length(lon),length(lat),1,dim(gpp)[3]))

## get time from a different file
time <- get_days_since( 1860, 1872, freq="months" )

cdf.write( gpp, "gpp",
           lon, lat,
           filnam = paste0( myhome, "data/trendy/v5/SDGVM/S2/SDGVM_S2_gpp_NICE.nc" ),
           nvars = 1,
           time = time,
           make.zdim = TRUE,
           z_dim = 1,
           make.tdim = TRUE,
           units_time = "days since 1970-01-01 00:00:00",
           units_var1 = "kg C m-2 s-1",
           glob_hist = "created by soilm_global/preproc_sdgvm.R"
)
