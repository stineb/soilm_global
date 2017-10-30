library(ncdf4)
library(abind)

print("getting gpp ...")
nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/SDGVM/S2/SDGVM_S2_gpp.nc" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time_units <- nc$dim$time$units
time <- nc$dim$time$vals
gpp <- ncvar_get( nc, varid="gpp" )

gpp_sub <- gpp[,,1465:1872]
time_sub <- nc$dim$time$vals[1465:1872]

## write reshaped GPP to file
cdf.write( gpp_sub, "gpp", 
           lon, lat,
           filnam = "/Users/benjaminstocker/data/trendy/v5/SDGVM/S2/SDGVM_S2_gpp_sub.nc",
           nvars = 1,
           time = time_sub,
           make.tdim = TRUE,
           units_time = time_units,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1"
)

nc_close( nc )
