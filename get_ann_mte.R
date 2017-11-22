source("~/.Rprofile")
library(dplyr)
library(ncdf4)
library(abind)

yearstart = 1982
yearend   = 2011

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
nsecsmonth <- ndaymonth * 60 * 60 * 24
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
vec_nsecsmonth <- rep( nsecsmonth, length(yearstart:yearend) )    ## This assumes that first time step is Jan 1982, but in file it says 21 Nov. 1981

fil <- paste0( myhome, "/data/gpp_mte/gpp_mte.nc")

print( paste("opening file: ", fil ) )
nc <- nc_open( fil )
gpp <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

gpp_permon <- sweep( gpp, 3, vec_nsecsmonth, "*" ) * 1e3  # somehow, it was to small by factor 1000 (probaly meant to be kg instead of g)

# ## This is a correct application of 'sweep' as demonstrated below...
# gpp_permon_test <- gpp * NA
# for (ilon in seq(length(lon))){
#   for (ilat in seq(length(lat))){
#     gpp_permon_test[ilon,ilat,] <- gpp[ilon,ilat,] * vec_nsecsmonth * 1e3
#   }
# }
# print(all.equal(gpp_permon, gpp_permon_test))


## write GPP in units of gC m-2 month-1 to file
cdf.write( gpp_permon, "gpp", 
           lon, lat,
           filnam = paste0(myhome, "/data/gpp_mte/gpp_mte_permon.nc"),
           nvars = 1,
           time = time,
           make.tdim = TRUE,
           units_time = "month",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1",
           glob_hist = "file created by utilities/get_ann_mte.nc based on data/gpp_mte/gpp_mte.nc, received from T. Keenan."
)

## Reshape array to have a separate dimension for months
gpp_resh <- array( gpp_permon, dim=c( length(lon), length(lat), nmonth, dim(gpp_permon)[3]/nmonth ) )

## Take annual GPP
gpp_ann <- apply( gpp_resh, c(1,2,4), FUN = sum )

## write annual GPP to file
outfilnam <- paste0(myhome, "/data/gpp_mte/gpp_mte_ann.nc")
cdf.write( gpp_ann, "gpp", 
           lon, lat,
           filnam = outfilnam,
           nvars = 1,
           time = yearstart:yearend,
           make.tdim = TRUE,
           units_time = "year",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1",
           glob_hist = "file created by utilities/get_ann_mte.nc based on data/gpp_mte/gpp_mte.nc, received from T. Keenan."
)



## ALTERNATIVELY, USE DATA DOWNLOADED FROM FLUXCOM WEBSITE:
yearstart = 1980
yearend   = 2013

idx <- 0
for (year in yearstart:yearend){

  idx <- idx + 1
  filn <- paste0( myhome, "/data/gpp_mte/GPP.RF.CRUNCEPv6.annual.", as.character(year), ".nc" )

  nc <- nc_open( filn )
  tmp <- ncvar_get( nc, varid="GPP" )
  tmp <- tmp * 365  # convert to totals, given in units per day
  if (idx>1) {
    gpp <- abind( gpp, tmp, along = 3 )
  } else {
    gpp <- tmp
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
  }
  nc_close(nc)

}

## write annual GPP to file
outfilnam <- paste0(myhome, "/data/gpp_mte/gpp_mte_fluxcom_ann.nc")
cdf.write( gpp, "gpp", 
           lon, lat,
           filnam = outfilnam,
           nvars = 1,
           time = yearstart:yearend,
           make.tdim = TRUE,
           units_time = "year",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1",
           glob_hist = "file created by utilities/get_ann_mte.nc based on GPP.RF.CRUNCEPv6.annual.<year>.nc, downloaded from ftp://ftp.bgc-jena.mpg.de/pub/outgoing/FluxCom/CarbonFluxes/RS+METEO/CRUNCEPv6/raw/annual/ (17.11.2017)."
)
