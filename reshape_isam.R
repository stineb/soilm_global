library(ncdf4)
library(abind)

print("getting gpp ...")
nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/ISAM/S2/ISAM_S2_gpp_sub_tmp.nc" )
lon <- nc$dim$LON$vals
lat <- nc$dim$LAT$vals
gpp <- ncvar_get( nc, varid="GPP" )

print("getting time from CLM file ...")
nc_time <- nc_open( "/Users/benjaminstocker/data/trendy/v5/CLM/S2/CLM4.5_S2_gpp_sub.nc" )
time_units <- nc_time$dim$time$units
time <- nc_time$dim$time$vals

gpp_resh <- array( gpp, dim=c(length(lon), length(lat), dim(gpp)[4]*12) )

## write reshaped GPP to file
cdf.write( gpp_resh, "gpp", 
           lon, lat,
           filnam = "/Users/benjaminstocker/data/trendy/v5/ISAM/S2/ISAM_S2_gpp_sub.nc",
           nvars = 1,
           time = time+31,
           make.tdim = TRUE,
           units_time = time_units,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1"
)

nc_close( nc )
nc_close( nc_time )
