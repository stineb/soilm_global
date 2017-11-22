source("~/.Rprofile")
library(dplyr)
library(ncdf4)
library(pracma)   # provides function 'detrend'

source("integrate_gridcell.R")

## Files for which variance is derived from 30 years data (1982-2011)
filnams_30y <- c( paste0( myhome, "/data/gpp_mte/gpp_mte"),
                  paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom")
                  )
time_30y <- 1982:2011

## Files for which variance is derived from 13 years data (2001-2013)
filnams_10y <- c( paste0( myhome, "/data/gpp_bess/gpp_bess"),
                  paste0( myhome, "/data/gpp_mte/gpp_mte"),
                  paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom"),
                  paste0( myhome, "data/gpp_modis/gpp_modis"),
                  paste0( myhome, "data/gpp_vpm/gpp_vpm")
                  )

## WARNING: WILL HAVE TO EXTEND THIS TO 2013 (AT LEAST) AS SOON AS P-MODEL SIMULATION IS DONE WITH LONGER TIME SERIES
time_10y <- 2001:2011

filnams <- list( l10y=filnams_10y, l30y=filnams_30y )

for (ilen in c("l30y", "l10y")){

  print( paste( "Getting variance for length", ilen ) )

  for ( basefil in filnams[[ ilen ]] ){

    if (ilen=="l30y"){
      ext <- ""
    } else if (ilen=="l10y"){
      ext <- "_10y"
    }

    if (!file.exists(paste0( myhome, basefil, "_relvar", ext, ".nc"))){

      print( paste("opening file: ", paste0(basefil, "_ann.nc") ) )
      nc <- nc_open( paste0(basefil, "_ann.nc") )
      gpp <- ncvar_get( nc, varid="gpp" )
      lon <- nc$dim$lon$vals
      lat <- nc$dim$lat$vals
      time <- nc$dim$time$vals
      nc_close(nc)

      dx <- abs(lon[2] - lon[1])
      dy <- abs(lat[2] - lat[1])

      print("found years:")
      print(range(time))

      ## reducing years
      if (ext==""){
        time <- time[which(time %in% time_30y)]
        gpp <- gpp[,,which(time %in% time_30y)]
      } else {
        time <- time[which(time %in% time_10y)]
        gpp <- gpp[,,which(time %in% time_10y)]
      }
      
      gpp_nice <- gpp

      print("integrating globally")
      ggpp <-  integrate_gridcell( gpp_nice, overwrite = FALSE ) * 1e-15

      ## write dataframe of global totals
      print("save global total timeseries")
      df <- data.frame( year=time, gpp=ggpp )
      write.csv( df, file=paste0( basefil, "_globaltotal", ext, ".csv" ), row.names=FALSE )
      save( df, file=paste0( basefil, "_globaltotal", ext, ".Rdata" ) )

      ## store nice and unified GPP
      print("save nice file")
      outfilnam <- paste0( basefil, "_nice", ext, ".nc")
      cdf.write( gpp_nice, "gpp", 
                 lon, lat,
                 filnam = outfilnam,
                 nvars = 1,
                 time = time,
                 make.tdim = TRUE,
                 units_time = "year",
                 long_name_var1 = "Gross primary productivity",
                 units_var1 = "gC m-2 year-1"
      )

      ## Detrend data
      print("detrend data")
      gpp_detr <- apply( gpp_nice, c(1,2), FUN = detrend )
      gpp_detr <- aperm( gpp_detr, c(2,3,1) )
      cdf.write( gpp_detr, "gpp", 
                 lon, lat,
                 filnam = paste0( basefil, "_detr", ext, ".nc"),
                 nvars = 1,
                 time = time,
                 make.tdim = TRUE,
                 units_time = "year",
                 long_name_var1 = "Gross primary productivity",
                 units_var1 = "gC m-2 year-1"
      )

      ## Get mean, variance, and relative variance
      print("get variance")
      gpp_var <- apply( gpp_detr, c(1,2), FUN = var )
      cdf.write( gpp_var, "gpp", 
                 lon, lat,
                 filnam = paste0( basefil, "_var", ext, ".nc"),
                 nvars = 1,
                 make.tdim = FALSE,
                 long_name_var1 = "Gross primary productivity",
                 units_var1 = "gC m-2 year-1"
      )

      print("get mean")
      gpp_mean <- apply( gpp_nice, c(1,2), FUN = mean )
      cdf.write( gpp_mean, "gpp", 
                 lon, lat,
                 filnam = paste0( basefil, "_mean", ext, ".nc"),
                 nvars = 1,
                 make.tdim = FALSE,
                 long_name_var1 = "Gross primary productivity",
                 units_var1 = "gC m-2 year-1"
      )

      print("get relative variance")
      gpp_relvar <- gpp_var / gpp_mean
      cdf.write( gpp_relvar, "gpp", 
                 lon, lat,
                 filnam = paste0( basefil, "_relvar", ext, ".nc"),
                 nvars = 1,
                 make.tdim = FALSE,
                 long_name_var1 = "Gross primary productivity",
                 units_var1 = "(unitless)"
      )

    }

  }

}


## Analyse agreement of variance when derived from 10 years and from 30 years
source( paste0( myhome, "sofun/utils_sofun/analysis_sofun/analyse_modobs.R" ) )
for ( basefil in filnams[[ "l30y" ]] ){

  print( paste("opening file: ", paste0(basefil, "_var.nc") ) )
  nc <- nc_open( paste0(basefil, "_var.nc") )
  var_30y <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)


  print( paste("opening file: ", paste0(basefil, "_var_10y.nc") ) )
  nc <- nc_open( paste0(basefil, "_var_10y.nc") )
  var_10y <- ncvar_get( nc, varid="gpp" )
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  time <- nc$dim$time$vals
  nc_close(nc)
  
  modobs <- analyse_modobs( c(var_10y), c(var_30y), plot.col =rgb(0,0,0,0.1) )
  
}



