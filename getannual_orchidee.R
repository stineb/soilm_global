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

gpp <- abind( gpp_1, gpp_2, gpp_3, gpp_4, along=3 )

# gpp_resh <- array( gpp, dim=c(length(lon), length(lat), 12, dim(gpp)[3]/12) )

# gpp_ann <- apply( gpp_resh, c(1,2,4), FUN = sum )

# ## write reshaped GPP to file
# cdf.write( gpp_ann, "gpp", 
#            lon, lat,
#            filnam = "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_ann.nc",
#            nvars = 1,
#            time = 1981:2016,
#            make.tdim = TRUE,
#            units_time = "years",
#            long_name_var1 = "Gross primary productivity",
#            units_var1 = "gC m-2 year-1"
# )

# cdf.write( gpp_ann[,,2:35], "gpp", 
#            lon, lat,
#            filnam = "/Users/benjaminstocker/data/trendy/v5/ORCHIDEE/S2/orchidee_S2_gpp_ann_sub.nc",
#            nvars = 1,
#            time = 1982:2015,
#            make.tdim = TRUE,
#            units_time = "years",
#            long_name_var1 = "Gross primary productivity",
#            units_var1 = "gC m-2 year-1"
# )