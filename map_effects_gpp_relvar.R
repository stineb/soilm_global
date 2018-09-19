library(ncdf4)
library(RColorBrewer)
library( ncdf4, quietly = TRUE )
library( fields, quietly = TRUE )
library( sp, quietly = TRUE )
library( maptools, quietly = TRUE )
library( dplyr, quietly = TRUE )  

source("../utilities/myboxplot.R")
source("plot_map.R")

## file name for figure
filn <- "fig/map_gpp_relvar_diff.pdf"

##------------------------------------------------------------------------
## GPP interannual relative variance change
##------------------------------------------------------------------------
fil_s0  <- "gpp_pmodel_s0_RELVAR.nc"
fil_s1  <- "gpp_pmodel_s1_RELVAR.nc"
fil_s1a <- "gpp_pmodel_s1a_RELVAR.nc"
fil_s1b <- "gpp_pmodel_s1b_RELVAR.nc"
fil_s1c <- "gpp_pmodel_s1c_RELVAR.nc"

dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

## S0
nc <- nc_open( paste0( dir, fil_s0 ) )
gpp_s0 <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

# ## S1
# nc <- nc_open( paste0( dir, fil_s1 ) )
# gpp_s1 <- ncvar_get( nc, varid="gpp" )
# nc_close(nc)
# 
# ## S1a
# nc <- nc_open( paste0( dir, fil_s1a ) )
# gpp_s1a <- ncvar_get( nc, varid="gpp" )
# nc_close(nc)

## S1b
nc <- nc_open( paste0( dir, fil_s1b ) )
gpp_s1b <- ncvar_get( nc, varid="gpp" )
nc_close(nc)

# ## S1c
# nc <- nc_open( paste0( dir, fil_s1c ) )
# gpp_s1c <- ncvar_get( nc, varid="gpp" )
# nc_close(nc)

## calculate amplification factor (field)
ampl <- gpp_s1b / gpp_s0

# ##-----------------------------------------------------
# ## Plot relative variance in S1
# ##-----------------------------------------------------
# # plot_map( gpp_s1b / gpp_s0, lev=c( 0, 4, 10 ),
# #       toplefttext=expression(paste("Amplification of GPP relative variance")),
# #       toprighttext=expression(paste("fraction")),
# #       maxval = 35, positive = FALSE, color = c( "royalblue4", "wheat", "tomato2", "tomato4" )
# #       )

# plot_map( gpp_s1b - gpp_s0, lev=c( -5, 15, 10 ),
#           toplefttext=expression(paste("Difference in GPP relative variance")),
#           toprighttext=expression(paste("unitless")),
#           maxval = 35, positive = FALSE, color = c( "royalblue4", "wheat", "tomato2", "tomato4" )
#         )

# plot_map( gpp_s1b - gpp_s0, lev=c( -7, 13, 10 ),
#           toplefttext=expression(paste("Difference in GPP relative variance")),
#           toprighttext=expression(paste("unitless")),
#           maxval = 35, positive = FALSE, color = c( "royalblue4", "wheat", "tomato2", "tomato4" )
#         )

# ## Density distribution of relative variance in s0 and s1b
# hist( gpp_s0, xlim=c(0,40), breaks=300, col=rgb(0,0,0,0.3) )
# hist( gpp_s1b, breaks=300, col=rgb(1,0,0,0.3), add=TRUE )


# hist(ampl, xlim=c(0,10), breaks=300)
# abline(v=1, col="red")

## Analyse distribution factor vs. aridity (mean annual AET/PET)
ncfiln <- "../data/greve/ep_over_p_cru_ncep.nc"
if (!file.exists(ncfiln)) {
  epop <- array( 1, dim=c(720,360) )
} else {
  nc <- nc_open( ncfiln )
  epop <- ncvar_get( nc, varid="EP_OVER_P_CRU_NCEP" )
}
alpha <- 1/epop

# hist( 1/epop, xlim=c(1,5), breaks = 3000 )

df_pixels <-  tibble( ampl = c(ampl), alpha = c(alpha) ) %>% 
		      filter( !is.na(ampl) & !is.na(alpha) ) %>% 
		      mutate( inbin = cut( alpha, breaks = c(0, 0.05, 0.2, 0.5, 0.7, 1.3, 3) ) )


