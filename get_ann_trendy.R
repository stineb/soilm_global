source("~/.Rprofile")
library(dplyr)

ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
nsecsmonth <- ndaymonth * 60 * 60 * 24
ndayyear <- sum(ndaymonth)
nmonth <- length(ndaymonth)
vec_nsecsmonth <- rep( nsecsmonth, length(1982:2015) )

if (!exists("filnams")) filnams <- read.csv("/Users/benjaminstocker/data/trendy/v5/trendy_s2_filnams_gpp.csv")

# for (idx in 8:nrow(filnams)){
for (idx in 8:8){

	if (filnams$orig[idx]!=""){

		modl <- as.character(filnams$modl[idx])
		fil <- as.character(filnams$sub[idx])
		basefil <- substr(fil, start=1, stop=nchar(fil)-3)

		print( paste("converting to annual total: ", fil ) )
		nc <- nc_open( paste0("/Users/benjaminstocker/data/trendy/v5/", modl, "/S2/", fil ) )
		gpp <- ncvar_get( nc, varid="gpp" )
		lon <- nc$dim[[ as.character(filnams$lonname[idx]) ]]$vals
		lat <- nc$dim[[ as.character(filnams$latname[idx]) ]]$vals
		time_units <- nc$dim[[ as.character(filnams$timename[idx]) ]]$units
		time <- nc$dim[[ as.character(filnams$timename[idx]) ]]$vals
		nc_close(nc)

		gpp_permon <- sweep( gpp, 3, vec_nsecsmonth, "*" )

		## write GPP in units of gC m-2 month-1 to file
		cdf.write( gpp_permon, "gpp", 
		           lon, lat,
		           filnam = paste0("/Users/benjaminstocker/data/trendy/v5/", modl, "/S2/", basefil, "_permon.nc"),
		           nvars = 1,
		           time = time,
		           make.tdim = TRUE,
		           units_time = time_units,
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 month-1"
		)

		## Reshape array to have a separate dimension for months
		gpp_resh <- array( gpp_permon, dim=c( length(lon), length(lat), nmonth, dim(gpp_permon)[3]/nmonth ) )

		## Take annual GPP
		gpp_ann <- apply( gpp_resh, c(1,2,4), FUN = sum )

		## write annual GPP to file
		outfilnam <- paste0("/Users/benjaminstocker/data/trendy/v5/", modl, "/S2/", basefil, "_ann.nc")
		cdf.write( gpp_ann, "gpp", 
		           lon, lat,
		           filnam = outfilnam,
		           nvars = 1,
		           time = 1982:2015,
		           make.tdim = TRUE,
		           units_time = "year",
		           long_name_var1 = "Gross primary productivity",
		           units_var1 = "gC m-2 year-1"
		)
		filnams$ann[idx] <- outfilnam

	}

}

write.csv( filnams, file="/Users/benjaminstocker/data/trendy/v5/trendy_s2_filnams_gpp.csv" )