source("~/.Rprofile")
library(dplyr)
library(pracma)

filnams <- read.csv("/Users/benjaminstocker/data/trendy/v5/trendy_s2_filnams_gpp.csv")
filnams$nice <- rep(NA, nrow(filnams))

# for (idx in 1:nrow(filnams)){
for (idx in 3:3){

	if (filnams$orig[idx]!=""){

		modl <- as.character(filnams$modl[idx])
		fil <- as.character(filnams$ann[idx])
		basefil <- substr(fil, start=1, stop=nchar(fil)-3)

		print( paste("opening file: ", fil ) )
		nc <- nc_open( fil )
		gpp <- ncvar_get( nc, varid="gpp" )
		lon <- nc$dim$LONGITUDE$vals
		lat <- nc$dim$LATITUDE$vals
		time <- nc$dim$TIME$vals
		nc_close(nc)

		dx <- abs(lon[2] - lon[1])
		dy <- abs(lat[2] - lat[1])

		if (modl=="SDGVM"){
			## SDGVM is seriously messed up
			lon <- lon - 180
			lat <- lat - 90
			gpp <- abs( gpp )
		}

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

	  print("integrating globally")
	  if (modl=="ISAM" || modl=="SDGVM"){
	  	## assuming that ISAM and SDGVM provide gridcell integrated values already, otherwise global totals do not make sense
	  	gpp_abs <- gpp
	  	gpp_nice <- sweep( gpp, c(1,2), arr_area, "/" )
	  } else {
		  gpp_abs <- sweep( gpp, c(1,2), arr_area, "*" )
		  gpp_nice <- gpp
	  }

	  ggpp <- apply( gpp_abs, c(3), FUN = sum, na.rm=TRUE )

	  if (modl=="CABLE" || modl=="CLM" || modl=="JSBACH" || modl=="LPJ-GUESS" || modl=="LPX-Bern" || modl=="ORCHIDEE" || modl=="VEGAS" || modl=="VISIT"){
	  	ggpp <- ggpp * 1e-12
	  	gpp_nice <- gpp_nice * 1e3
	  } else if (modl=="CLASS-CTEM"){
	  	ggpp <- ggpp * NA
	  	gpp_nice <- gpp_nice * NA
	  } else if (modl=="ISAM"){
	  	ggpp <- ggpp * 1e-9
	  	gpp_nice <- gpp_nice * 1e6
	  } else if (modl=="JULES"){
	  	ggpp <- ggpp * NA
	  	gpp_nice <- gpp_nice * NA
	  } else if (modl=="SDGVM"){
	  	ggpp <- ggpp * 1e-2
	  	gpp_nice <- gpp_nice * 1e13
	  }

	  ## write dataframe of global totals
	  print("save global total timeseries")
	  df <- data.frame( year=time, gpp=ggpp )
	  write.csv( df, file=paste0( "data/", modl, "_globaltotal.csv" ), row.names=FALSE )
	  save( df, file=paste0( "data/", modl, "_globaltotal.Rdata" ) )

	  print(paste("Model, range of global GPP:", modl, range(ggpp)))

	  ## store nice and unified GPP
	  print("save nice file")
		outfilnam <- paste0( basefil, "_nice.nc" )
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
		filnams$nice[idx] <- outfilnam

		## Detrend data
	  print("detrend data")
		gpp_detr <- apply( gpp_nice, c(1,2), FUN = detrend )
		gpp_detr <- aperm( gpp_detr, c(2,3,1) )
		cdf.write( gpp_detr, "gpp", 
		           lon, lat,
		           filnam = paste0( basefil, "_detr.nc" ),
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
		           filnam = paste0( basefil, "_var.nc" ),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

	  print("get mean")
		gpp_mean <- apply( gpp_nice, c(1,2), FUN = mean )
		cdf.write( gpp_mean, "gpp", 
		           lon, lat,
		           filnam = paste0( basefil, "_mean.nc" ),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

	  print("get relative variance")
		gpp_relvar <- gpp_var / gpp_mean
		cdf.write( gpp_relvar, "gpp", 
		           lon, lat,
		           filnam = paste0( basefil, "_relvar.nc" ),
		           nvars = 1,
		           make.tdim = FALSE,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)

	}

}

