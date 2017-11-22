source("~/.Rprofile")
library(dplyr)
library(ncdf4)
library(pracma)   # provides function 'detrend'

source("integrate_gridcell.R")

simnames <- c( "s0_fapar3g_global", "s1_fapar3g_global" )
dir <- paste0( myhome, "/data/pmodel_fortran_output/" )

time_30y <- 1982:2011

## WARNING: NOW ONLY TO 2011 - SHOULD BE TO 2013
time_10y <- 2001:2011

for (ilen in c("l30y", "l10y")){

	for (ifil in simnames){

    if (ilen=="l30y"){
      ext <- ""
    } else if (ilen=="l10y"){
      ext <- "_10y"
    }

		print( paste("opening file: ", ifil ) )
		nc <- nc_open( paste0( dir, ifil, ".a.gpp.nc" ) )
		gpp <- ncvar_get( nc, varid="gpp" )
		lon <- nc$dim$lon$vals
		lat <- nc$dim$lat$vals
		time <- nc$dim$time$vals
		nc_close(nc)

		## WARNING: HACK
		time <- 1982:2011

		dx <- abs(lon[2] - lon[1])
		dy <- abs(lat[2] - lat[1])

    ## reducing years
    if (ext==""){
      time <- time[which(time %in% time_30y)]
      gpp  <- gpp[,,which(time %in% time_30y)]
      time_out <- time_30y
    } else {
      time <- time[which(time %in% time_10y)]
      gpp  <- gpp[,,which(time %in% time_10y)]
      time_out <- time_10y
    }

		gpp_nice <- gpp

		print("integrating globally")
		ggpp <- integrate_gridcell( gpp_nice )
		ggpp <- ggpp * 1e-15

		## write dataframe of global totals
		print("save global total timeseries")
		df <- data.frame( year=time_out, gpp=ggpp )
		write.csv( df, file=paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_globaltotal", ifil , ext, ".csv" ), row.names=FALSE )
		save( df, file=paste0( "data/pmodel_gpp_globaltotal", ifil , ext, ".Rdata" ) )

		print(paste("Model, range of global GPP: P-model, simulation  ", ifil, range(ggpp)))

		## store nice and unified GPP
		print("save nice file")
		outfilnam <- paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_nice_", ifil , ext, ".nc")
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
		           filnam = paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_detr_", ifil , ext, ".nc"),
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
		           filnam = paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_var_", ifil , ext, ".nc"),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

		print("get mean")
		gpp_mean <- apply( gpp_nice, c(1,2), FUN = mean )
		cdf.write( gpp_mean, "gpp", 
		           lon, lat,
		           filnam = paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_mean_", ifil , ext, ".nc"),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

		print("get relative variance")
		gpp_relvar <- gpp_var / gpp_mean
		cdf.write( gpp_relvar, "gpp", 
		           lon, lat,
		           filnam = paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_relvar_", ifil , ext, ".nc"),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

	}

}

## Analyse agreement of variance when derived from 10 years and from 30 years
source( paste0( myhome, "sofun/utils_sofun/analysis_sofun/analyse_modobs.R" ) )

for (ifil in simnames){

	print( paste("opening file: ", ifil ) )
	nc <- nc_open( paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_var_", ifil, ".nc") )
	var_30y <- ncvar_get( nc, varid="gpp" )
	lon <- nc$dim$lon$vals
	lat <- nc$dim$lat$vals
	time <- nc$dim$time$vals
	nc_close(nc)

	print( paste("opening file: ", ifil ) )
	nc <- nc_open( paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_var_", ifil , "_10y.nc") )
	var_10y <- ncvar_get( nc, varid="gpp" )
	lon <- nc$dim$lon$vals
	lat <- nc$dim$lat$vals
	time <- nc$dim$time$vals
	nc_close(nc)

	modobs <- analyse_modobs( c(var_10y), c(var_30y), plot.col =rgb(0,0,0,0.1) )

}


