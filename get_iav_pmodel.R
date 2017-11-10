source("~/.Rprofile")
library(dplyr)
library(ncdf4)
library(pracma)   # provides function 'detrend'

simnames <- c( "s0_fapar3g_global", "s1_fapar3g_global" )
dir <- "/Users/benjaminstocker/sofun/trunk/output_nc/"

for (ifil in simnames){

	print( paste("opening file: ", ifil ) )
	nc <- nc_open( paste0( dir, ifil, ".a.gpp.nc" ) )
	gpp <- ncvar_get( nc, varid="gpp" )
	lon <- nc$dim$lon$vals
	lat <- nc$dim$lat$vals
	time <- nc$dim$time$vals
	nc_close(nc)

	dx <- abs(lon[2] - lon[1])
	dy <- abs(lat[2] - lat[1])

	print("getting area array")
	arr_area <- gpp[,,1]
	arr_area[] <- NA
	for (ilon in seq(dim(gpp)[1])){
	  for (ilat in seq(dim(gpp)[2])){
	    if (!is.na(gpp[ilon,ilat,])){
	      arr_area[ilon,ilat] <- area( lat[ilat], dx=dx, dy=dy )
	    }
	  }
	}

	gpp_nice <- gpp

	print("integrating globally")
	gpp_abs <- sweep( gpp_nice, c(1,2), arr_area, "*" )
	ggpp <- apply( gpp_abs, c(3), FUN = sum, na.rm=TRUE )
	ggpp <- ggpp * 1e-15

	## write dataframe of global totals
	print("save global total timeseries")
	df <- data.frame( year=1982:2011, gpp=ggpp )
	write.csv( df, file=paste0( "data/pmodel_gpp_globaltotal", ifil ,".csv" ), row.names=FALSE )
	save( df, file=paste0( "data/pmodel_gpp_globaltotal", ifil ,".Rdata" ) )

	print(paste("Model, range of global GPP: P-model, simulation  ", ifil, range(ggpp)))

	## store nice and unified GPP
	print("save nice file")
	outfilnam <- paste0( "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_nice_", ifil ,".nc")
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
	           filnam = paste0( "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_detr_", ifil ,".nc"),
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
	           filnam = paste0( "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_var_", ifil ,".nc"),
	           nvars = 1,
	           make.tdim = FALSE,
	           long_name_var1 = "Gross primary productivity",
	           units_var1 = "gC m-2 year-1"
	)

	print("get mean")
	gpp_mean <- apply( gpp_nice, c(1,2), FUN = mean )
	cdf.write( gpp_mean, "gpp", 
	           lon, lat,
	           filnam = paste0( "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_mean_", ifil ,".nc"),
	           nvars = 1,
	           make.tdim = FALSE,
	           long_name_var1 = "Gross primary productivity",
	           units_var1 = "gC m-2 year-1"
	)

	print("get relative variance")
	gpp_relvar <- gpp_var / gpp_mean
	cdf.write( gpp_relvar, "gpp", 
	           lon, lat,
	           filnam = paste0( "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_relvar_", ifil ,".nc"),
	           nvars = 1,
	           make.tdim = FALSE,
	           long_name_var1 = "Gross primary productivity",
	           units_var1 = "gC m-2 year-1"
	)

}
