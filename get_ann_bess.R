source("~/.Rprofile")
library(dplyr)
library(ncdf4)

yearstart = 2001
yearend   = 2015

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
vec_ndaymonth <- rep( ndaymonth, length(yearstart:yearend) )    ## This assumes that first time step is Jan 1982, but in file it says 21 Nov. 1981

fil <- "/Users/benjaminstocker/data/gpp_bess/BESS_GPP_Monthly_05D.nc"

print( paste("opening file: ", fil ) )
nc <- nc_open( fil )
gpp <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

gpp_permon <- sweep( gpp, 3, vec_ndaymonth, "*" ) # somehow, it was to small by factor 1000 (probaly meant to be kg instead of g) 

## write GPP in units of gC m-2 month-1 to file
cdf.write( gpp_permon, "gpp", 
           lon, lat,
           filnam = "/Users/benjaminstocker/data/gpp_bess/gpp_bess_permon.nc",
           nvars = 1,
           time = time,
           make.tdim = TRUE,
           units_time = "days since 1970-01-01",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1",
           glob_hist = "created using soilm_global/get_ann_bess.R based on original file BESS_GPP_Monthly_05D.mat, downloaded from ftp://realtimedata.snu.ac.kr"
)

## Reshape array to have a separate dimension for months
gpp_resh <- array( gpp_permon, dim=c( length(lon), length(lat), nmonth, dim(gpp_permon)[3]/nmonth ) )

## Take annual GPP
gpp_ann <- apply( gpp_resh, c(1,2,4), FUN = sum )

## write annual GPP to file
outfilnam <- "/Users/benjaminstocker/data/gpp_bess/gpp_bess_ann.nc"
cdf.write( gpp_ann, "gpp", 
           lon, lat,
           filnam = outfilnam,
           nvars = 1,
           time = yearstart:yearend,
           make.tdim = TRUE,
           units_time = "year",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1",
           glob_hist = "created using soilm_global/get_ann_bess.R based on original file BESS_GPP_Monthly_05D.mat, downloaded from ftp://realtimedata.snu.ac.kr"
)
