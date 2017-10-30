library(ncdf4)
library(abind)

print("getting gpp ...")
nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_1981_1990.nc" )
lon <- nc$dim$longitude$vals
lat <- nc$dim$latitude$vals
gpp_1 <- ncvar_get( nc, varid="gpp" )
nc_close( nc )

nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_1991_2000.nc" )
gpp_2 <- ncvar_get( nc, varid="gpp" )
nc_close( nc )

nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_2001_2014.nc" )
gpp_3 <- ncvar_get( nc, varid="gpp" )
nc_close( nc )

nc <- nc_open( "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_2015_2016.nc" )
gpp_4 <- ncvar_get( nc, varid="gpp" )
nc_close( nc )

print("getting time from CLM file ...")
nc_time <- nc_open( "/Users/benjaminstocker/data/trendy/v5/CLM/S2/CLM4.5_S2_gpp_sub.nc" )
time_units <- nc_time$dim$time$units
time <- nc_time$dim$time$vals
# calendar <- nc_time$dim$time$calendar

gpp <- abind( gpp_1, gpp_2, gpp_3, gpp_4, along=3 )

gpp_sub <- gpp[,,13:(dim(gpp)[3]-12)]

cdf.write( gpp_sub, "gpp", 
           lon, lat,
           filnam = "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_sub.nc",
           nvars = 1,
           time = time+30,
           make.tdim = TRUE,
           units_time = time_units,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1"
)

nc_close( nc )
nc_close( nc_time )