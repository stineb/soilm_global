library(ncdf4)
library(RColorBrewer)
source("plot_map.R")

source("~/.Rprofile")

##------------------------------------------------------------------------
## GPP interannual (absolute) variance change
##------------------------------------------------------------------------
fil_s0 <- "gpp_pmodel_s0_VAR.nc"
fil_s1b <- "gpp_pmodel_s1b_VAR.nc"

dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

## S0
nc <- nc_open( paste0( dir, fil_s0 ) )
gpp_s0 <- ncvar_get( nc, varid="gpp" )
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
time <- nc$dim$time$vals
nc_close(nc)

## S1
nc <- nc_open( paste0( dir, fil_s1b ) )
gpp_s1b <- ncvar_get( nc, varid="gpp" )
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

# ##-----------------------------------------------------
# ## Amplification in absolute variance
# ##-----------------------------------------------------
# color <- c( "royalblue4", "wheat", "tomato2", "tomato4" )
# lev <- c( 0, 4, 10 )
# maxval = 35
# minval = NA
# # par( mar=c(4,3,3,1),xaxs="i", yaxs="i",las=1, mgp=c(3,1,0))
# out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=maxval, minval=minval )
# plot_map( gpp_s1b/gpp_s0, lev=lev, maxval=maxval, color = out.mycolorbar$colors )

# plot_map( gpp_s1b - gpp_s0, lev=seq(-5000,13000,2000), minval=-169000, maxval=83000, toplefttext=expression(paste("")), toprighttext=expression(paste("g C m"^-2, " yr"^-1 ) ), color=c( "royalblue4", "wheat", "tomato2", "tomato4" )  )

# hist( gpp_s0, xlim=c(0,50000), breaks=300, col=rgb(0,0,0,0.3) )
# hist( gpp_s1b, breaks=300, col=rgb(1,0,0,0.3), add=TRUE )

# ## Analyse distribution factor vs. aridity (mean annual AET/PET)
# ncfiln <- "../data/greve/ep_over_p_cru_ncep.nc"
# if (!file.exists(ncfiln)) {
#   epop <- array( 1, dim=c(720,360) )
# } else {
#   nc <- nc_open( ncfiln )
#   epop <- ncvar_get( nc, varid="EP_OVER_P_CRU_NCEP" )
# }
# alpha <- 1/epop

# diff <- gpp_s1b - gpp_s0
# df <- tibble( diff = c(diff), alpha = c(alpha) ) %>% 
#       filter( !is.na(diff) & !is.na(alpha) ) %>% 
#       mutate( inbin = cut( alpha, breaks = c(0, 0.05, 0.2, 0.5, 0.7, 1.3, 3) ) )

# par(las=1)
# myboxplot( log(diff) ~ inbin, data = df, outline=FALSE,
#            col=colorRampPalette( rev( c( "royalblue4", "wheat", "tomato2", "tomato4" ) ) )( 6 ),
#            ylab="Difference", xlab="Aridity index")
# abline(h=1, lty=3)      
