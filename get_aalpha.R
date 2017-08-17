library(ncdf4)
library(dplyr)

## get AET
if (!exists("aet")){
  print("getting aet ...")
  nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_AET_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  aet <- ncvar_get( nc, varid="AET" )
  lon <- nc$dim$longitude$vals
  lat <- nc$dim$latitude$vals
  nc_close( nc )
}

## get PET
if (!exists("pet")){
  print("getting pet ...")
  nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_PET_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  pet <- ncvar_get( nc, varid="PET" )
  nc_close( nc )
}

## calculate monthly alpha = monthly AET / monthly PET
if (!exists("malpha")) malpha <- aet / pet
malpha[ which(is.infinite(malpha)) ] <- NA

## get mean AET/PET as an average across annual values
if (!exists("aalpha")) aalpha <- apply( malpha, c(1,2), FUN = function(x) min( 1.0, mean( x, na.rm=TRUE ) ) )

## write to file
cdf.write(  aalpha, "alpha", 
            lon, lat, 
            filnam = "data/alpha_Pmod_SPLASH.nc", 
            long_name_var1 = "annual mean AET/PET, based on monthly values", 
            units_var1 = "fraction" 
            )
