source("~/.Rprofile")
library(dplyr)
library(ncdf4)
library(pracma)   # provides function 'detrend'

source("integrate_gridcell.R")

fil <- paste0( myhome, "/data/gpp_mte/gpp_mte_ann.nc")

print( paste("opening file: ", fil ) )
nc <- nc_open( fil )
gpp <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$LONGITUDE$vals
lat <- nc$dim$LATITUDE$vals
time <- nc$dim$TIME$vals
nc_close(nc)

dx <- abs(lon[2] - lon[1])
dy <- abs(lat[2] - lat[1])

gpp_nice <- gpp * 1e3

print("integrating globally")
ggpp <- integrate_gridcell( gpp_nice )
ggpp <- ggpp * 1e-15

## write dataframe of global totals
print("save global total timeseries")
df <- data.frame( year=time, gpp=ggpp )
write.csv( df, file=paste0( "data/mte_globaltotal.csv" ), row.names=FALSE )
save( df, file=paste0( "data/mte_globaltotal.Rdata" ) )

print(paste("Model, range of global GPP: MTE  ", range(ggpp)))

## store nice and unified GPP
print("save nice file")
outfilnam <- paste0( myhome, "/data/gpp_mte/MTE_nice.nc")
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
           filnam = paste0( myhome, "/data/gpp_mte/MTE_detr.nc"),
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
           filnam = paste0( myhome, "/data/gpp_mte/MTE_var.nc"),
           nvars = 1,
           make.tdim = FALSE,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1"
)

print("get mean")
gpp_mean <- apply( gpp_nice, c(1,2), FUN = mean )
cdf.write( gpp_mean, "gpp", 
           lon, lat,
           filnam = paste0( myhome, "/data/gpp_mte/MTE_mean.nc"),
           nvars = 1,
           make.tdim = FALSE,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1"
)

print("get relative variance")
gpp_relvar <- gpp_var / gpp_mean
cdf.write( gpp_relvar, "gpp", 
           lon, lat,
           filnam = paste0( myhome, "/data/gpp_mte/MTE_relvar.nc"),
           nvars = 1,
           make.tdim = FALSE,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1"
)
