library(ncdf4)
library(RColorBrewer)
source("plot_map.R")
source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP interannual (absolute) variance change
##------------------------------------------------------------------------
fil_s0 <- "gpp_pmodel_s0_VAR.nc"
fil_s1 <- "gpp_pmodel_s1b_VAR.nc"

dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

## S0
nc <- nc_open( paste0( dir, fil_s0 ) )
gpp_s0 <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

## S1
nc <- nc_open( paste0( dir, fil_s1 ) )
gpp_s1 <- ncvar_get( nc, varid="gpp" )
nc_close(nc)

# ##-----------------------------------------------------
# ## Plot absolute variance in S1
# ##-----------------------------------------------------
# plot_map( gpp_s1*1e-3, lev=c( 0, 40, 10 ), maxval = 200 )
# 
# ##-----------------------------------------------------
# ## Difference in absolute variance
# ##-----------------------------------------------------
# plot_map( (gpp_s1 - gpp_s0)*1e-3, lev=c( -10, 10, 10 ), positive = FALSE, maxval=60, minval=-60 )

##-----------------------------------------------------
## Amplification in absolute variance
##-----------------------------------------------------
color <- c( "royalblue4", "wheat", "tomato2", "tomato4" )
lev <- c( 0, 4, 10 )
maxval = 35
minval = NA
# par( mar=c(4,3,3,1),xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=maxval, minval=minval )

plot_map( gpp_s1/gpp_s0, lev=lev, maxval=maxval, color = out.mycolorbar$colors )

