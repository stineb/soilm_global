library(ncdf4)
library(dplyr)
library(fields)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste0( myhome, "/utilities/mean.bymonth.R" ) )

## get annual mean alpha (AET/PET)
if (!exists("aalpha")){
  print("getting alpha ...")
  nc <- nc_open( "data/alpha_Pmod_SPLASH.nc" )
  aalpha <- ncvar_get( nc, varid="alpha" )
  lon <- nc$dim$longitude$vals
  lat <- nc$dim$latitude$vals
  nc_close( nc )
}

## get GPP
if (!exists("gpp")){
  print("getting gpp ...")
  nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_GPP_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  gpp <- ncvar_get( nc, varid="GPP" )
  nc_close( nc )
}

## get soil moisture
if (!exists("soilm")){
  print("getting soilm ...")
  nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_SoilM_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  soilm <- ncvar_get( nc, varid="SoilM" )
  soilm <- soilm / 150   ## convert into a fraction of water holding capacity
  lon <- nc$dim$longitude$vals
  lat <- nc$dim$latitude$vals
  time <- nc$dim$time$vals
  units_time <- nc$dim$time$units
  nc_close( nc )
}

## get area matrix
if (!file.exists("data/area.nc")){

  arr_area <- soilm[,,1]
  arr_area[] <- NA
  for (ilon in seq(dim(aalpha)[1])){
    for (ilat in seq(dim(aalpha)[2])){
      if (!is.na(aalpha[ilon,ilat])){
        arr_area[ilon,ilat] <- area( lat[ilat], dx=0.5, dy=0.5 )
      }
    }
  }

  ## write area to file
  cdf.write( arr_area, "area", 
             lon, lat,
             long_name_var1 = "gridcell surface area",
             units_var1 = "m2",
             filnam = "data/area.nc" 
             )

} else {

  print("getting area ...")
  nc <- nc_open( "data/area.nc" )
  arr_area <- ncvar_get( nc, varid="area" )
  nc_close( nc )  

}

## get flue0 : maximum fLUE reduction at low soil moisture
print("getting flue0 ...")
calc_flue0 <- function( alpha ){ 
  ## equation from correct_bysoilm.R
  y <- 0.1 + 0.92 * alpha
  y <- min( y, 1.0 )
  y <- max( y, 0.0 )
  return( y )
}
flue0 <- apply( aalpha, c(1,2), FUN = calc_flue0 )


## get beta : parameter defining quadratic correction function
print("getting beta ...")
calc_beta <- function( flue0 ){
  x0 <- 0.125
  x1 <- 0.75
  y <- (flue0 - 1.0) / (x0 - x1)^2
  return( y )
}
beta <- apply( flue0, c(1,2), FUN = calc_beta )


## get time varying fLUE as a function of monthly mean soil moisture
print("getting flue ...")
calc_flue <- function( soilm, beta ){
  if ( is.na(soilm) ){
    y <- NA
  } else {
    x1 <- 0.75
    if (soilm>x1) {
      y <- 1
    } else {
      y <- beta * ( soilm - x1 )^2 + 1.0
    }
  }
  return( y )
}

monthly2annual <- function( vec ){
  if (length(vec)%%12 !=0){
    print("cant do it. returning NA")
    return(NA)
  } else {
    nyears <- length(vec) / 12
    out <- c()
    for (iyr in seq(nyears)){
      out <- c( out, sum( ((iyr-1)*12+1):((iyr-1)*12+12) ) )
    }
    return(out)
  }
}

## calculate fLUE for all cells
flue <- soilm
flue[] <- NA
for (ilon in seq(dim(aalpha)[1])){
  for (ilat in seq(dim(aalpha)[2])){
    if (!is.na(aalpha[ilon,ilat])){
      flue[ilon,ilat,] <- sapply( soilm[ilon,ilat,], FUN = function(x) calc_flue( x, beta[ilon,ilat] ) )
    }
  }
}

## get global total unlimited GPP over time
gpp_abs <- sweep( gpp, 1, arr_area, "*", check.margin=FALSE )
ggpp <- apply( gpp_abs, c(3), FUN=sum, na.rm=TRUE )

# ## Ok, 'sweep' works
# gpp_abs_test <- gpp
# gpp_abs_test[] <- NA
# for (ilon in seq(dim(aalpha)[1])){
#   for (ilat in seq(dim(aalpha)[2])){
#     if (!is.na(aalpha[ilon,ilat])){
#       gpp_abs_test[ilon,ilat,] <- gpp[ilon,ilat,] * arr_area[ilon,ilat]
#     }
#   }
# }
# all( gpp_abs==gpp_abs_test )

## get global total water-limited GPP over time
gpp_lim <- gpp * flue
gpp_lim_abs <- sweep( gpp_lim, 1, arr_area, "*", check.margin=FALSE )
ggpp_lim <- apply( gpp_lim_abs, c(3), FUN=sum, na.rm=TRUE )

## get global gpp of gridcells with AET/PET < 0.8
gpp_dry <- gpp
gpp_lim_dry <- gpp_lim
for (ilon in seq(dim(gpp_dry)[1])){
  for (ilat in seq(dim(gpp_dry)[2])){
    if (!is.na(aalpha[ilon,ilat])){
      if (aalpha[ilon,ilat]>0.75){
        gpp_dry[ilon,ilat,] <- 0
        gpp_lim_dry[ilon,ilat,] <- 0
      }
    }
  }
}
gpp_dry_abs <- sweep( gpp_dry, 1, arr_area, "*", check.margin=FALSE )
ggpp_dry <- apply( gpp_dry_abs, c(3), FUN=sum, na.rm=TRUE )

gpp_lim_dry_abs <- sweep( gpp_lim_dry, 1, arr_area, "*", check.margin=FALSE )
ggpp_lim_dry <- apply( gpp_lim_dry_abs, c(3), FUN=sum, na.rm=TRUE )


## write to file, monthly GPP
cdf.write( gpp_lim, "gpp_lim", 
           lon, lat,
           filnam = "data/gpp_lim.nc",
           nvars = 3,
           var2 = gpp, varnam2 = "gpp",
           var3 = flue, varnam3 = "flue",
           time = time,
           make.tdim = TRUE,
           units_time = units_time,
           long_name_var1 = "Gross primary productivity, limited by soil moisture",
           long_name_var2 = "Gross primary productivity, unlimited by soil moisture",
           units_var1 = "gC m-2 month-1",
           units_var2 = "gC m-2 month-1"
           )


## REASHAPE write to file, monthly GPP
gpp_resh <- array( gpp, dim=c(length(lon), length(lat), 12, dim(gpp)[3]/12) )
flue_resh <- array( flue, dim=c(length(lon), length(lat), 12, dim(flue)[3]/12) )
gpp_lim_resh <- array( gpp_lim, dim=c(length(lon), length(lat), 12, dim(gpp_lim)[3]/12) )

## write reshaped GPP to file
cdf.write( gpp_resh, "gpp", 
           lon, lat,
           filnam = "data/gpp_RESH.nc",
           nvars = 1,
           time = 1982:2011,
           make.tdim = TRUE,
           z_dim = 1:12,
           make.zdim = TRUE,
           units_time = units_time,
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 month-1"
)

cdf.write( gpp_lim_resh, "gpp_lim", 
           lon, lat,
           filnam = "data/gpp_lim_RESH.nc",
           nvars = 3,
           var2 = gpp_resh, varnam2 = "gpp",
           var3 = flue_resh, varnam3 = "flue",
           time = 1982:2011,
           make.tdim = TRUE,
           z_dim = 1:12,
           make.zdim = TRUE,
           units_time = units_time,
           long_name_var1 = "Gross primary productivity, limited by soil moisture",
           long_name_var2 = "Gross primary productivity, unlimited by soil moisture",
           units_var1 = "gC m-2 month-1",
           units_var2 = "gC m-2 month-1"
           )

## get mean difference due to soil moisture limitation (absolute)
diff <- gpp_lim - gpp
diff_meanseason <- apply( diff, c(1,2), FUN=mean.bymonth, na.rm=TRUE )
diff_mean <- apply( diff_meanseason, c(2,3), FUN=sum, na.rm=TRUE )
diff_mean[ which(is.na(aalpha))] <- NA

# ##------------------------
# ## absolute GPP loss, map
# ##------------------------
#   library(fields)
#   library(sp)
#   library(maptools)

#   magn <- 4
#   ncols <- 2
#   nrows <- 1
#   widths <- rep(1.6*magn,ncols)
#   widths[2] <- 0.2*widths[1]
#   heights <- rep(magn,nrows)
#   order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

#   ylim <- c(-60,85)
#   lat.labels <- seq(-90, 90, 30)
#   lat.short  <- seq(-90, 90, 10)
#   lon.labels <- seq(-180, 180, 60)
#   lon.short  <- seq(-180, 180, 10)

#   a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
#   b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

#   pdf( "fig/map_gpploss_abs_byflue.pdf", width=sum(widths), height=sum(heights) )

#     panel <- layout(
#               order,
#               widths=widths,
#               heights=heights,
#               TRUE
#               )
#     # layout.show( panel )

#     ## Color key
#     # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
#     color <- c("wheat3","tomato2","tomato4")
#     lev <- c( 0, 500, 10 )
#     out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=1600 )

#     par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
#     image(
#             seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
#             -diff_mean,
#             ylim=c(-60,85), 
#             # zlim=range(lev), 
#             yaxt="n", xaxt="n",
#             col=out.mycolorbar$colors, breaks=out.mycolorbar$margins,
#             xlab="", ylab=""
#             )
#     map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )

#     axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
#     axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )

#     axis( 4, at=lat.labels, lab=F, lwd=1.5 )
#     axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )

#     axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
#     axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )

#     axis( 3, at=lon.labels, lab=F, lwd=1.5 )
#     axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )

#     ## Color key
#     par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
#     out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1600 )

#   dev.off()  


## construct dataframes with monthly and annual values
mdf <-  data.frame( 
                  date=as.Date( time-2, origin='0001-01-01' ),
                  gpp=ggpp * 1e-15,
                  gpp_lim=ggpp_lim * 1e-15,
                  gpp_dry=ggpp_dry * 1e-15,
                  gpp_lim_dry=ggpp_lim_dry * 1e-15
                  ) %>%
        mutate( moy = as.numeric( format( date, format="%m" ) ),
                year = as.numeric( format( date, format="%Y" ) )
               )

adf <- mdf %>%  group_by( year ) %>% 
                summarise( gpp=sum(gpp), gpp_lim=sum(gpp_lim), gpp_dry=sum(gpp_dry), gpp_lim_dry=sum(gpp_lim_dry) ) %>% 
                mutate( gpp_inc=gpp/gpp[1], gpp_inc_lim=gpp_lim/gpp_lim[1], gpp120=gpp/gpp[1]*120, gpp_lim120=gpp_lim/gpp_lim[1]*120 )

##------------------------
## get linear trends and add residuals to data frame
##------------------------
  lm_gpp     <- lm( gpp ~ year, data=adf )
  lm_gpp_lim <- lm( gpp_lim ~ year, data=adf )
  adf$res_gpp <- lm_gpp$residuals
  adf$res_gpp_lim <- lm_gpp_lim$residuals

  lm_gpp120     <- lm( gpp120 ~ year, data=adf )
  lm_gpp_lim120 <- lm( gpp_lim120 ~ year, data=adf )
  adf$res_gpp120 <- lm_gpp120$residuals
  adf$res_gpp_lim120 <- lm_gpp_lim120$residuals

  lm_gpp_inc     <- lm( gpp_inc ~ year, data=adf )
  lm_gpp_inc_lim <- lm( gpp_inc_lim ~ year, data=adf )
  adf$res_gpp_inc <- lm_gpp_inc$residuals
  adf$res_gpp_inc_lim <- lm_gpp_inc_lim$residuals

  lm_gpp_dry     <- lm( gpp_dry ~ year, data=adf )
  lm_gpp_lim_dry <- lm( gpp_lim_dry ~ year, data=adf )
  adf$res_gpp_dry <- lm_gpp_dry$residuals
  adf$res_gpp_lim_dry <- lm_gpp_lim_dry$residuals


##------------------------
## plot time series, absolute
##------------------------
# pdf( "fig/gpp_global_abs_lim_unlim.pdf", width=6, height=5 )
  par(las=1, mar=c(4,4,1,1))
  with( adf, plot(year, gpp, type="l", ylim=c(130,180), ylab="GPP (PgC/yr)", lwd=1.5 ) )
  abline( lm_gpp, col="black", lty=2 )
  cf <- coef(lm_gpp) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-1, adj=0.1 )

  with( adf, lines(year, gpp_lim, col="red", lwd=1.5 ) )
  abline( lm_gpp_lim, col="red", lty=2 )
  cf <- coef(lm_gpp_lim) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-2, adj=0.1, col='red' )
# dev.off()

##------------------------
## plot time series, 120
##------------------------
  par(las=1)
  with( adf, plot(year, gpp120, type="l", ylim=c(115,140), ylab="GPP (PgC/yr)" ) )
  abline( lm_gpp120, col="black", lty=2 )
  cf <- coef(lm_gpp120) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-1, adj=0.1 )

  with( adf, lines(year, gpp_lim120, col="red" ) )
  abline( lm_gpp_lim120, col="red", lty=2 )
  cf <- coef(lm_gpp_lim120) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-2, adj=0.1, col='red' )  


##------------------------
## plot time series, relative increase
##------------------------
  par(las=1)
  with( adf, plot(year, gpp_inc, type="l", ylim=c(0.95,1.2), ylab="change in GPP (frac., rel. to 1982)" ) )
  abline( lm_gpp_inc, col="black", lty=2 )
  cf <- coef(lm_gpp_inc) %>% round( 4 )
  eq <- paste0( "slope = ", abs(cf[2]), " yr-1 " )
  mtext( eq, line=-1, adj=0.1 )

  with( adf, lines(year, gpp_inc_lim, col="red" ) )
  abline( lm_gpp_inc_lim, col="red", lty=2 )
  cf <- coef(lm_gpp_inc_lim) %>% round( 4 )
  eq <- paste0( "slope = ", abs(cf[2]), " yr-1 " )
  mtext( eq, line=-2, adj=0.1, col='red' )


##------------------------
## plot flue effect
##------------------------
  with( adf, plot( year, gpp_lim / gpp, type="l" ) )


##------------------------
## plot residuals from linear trend
##------------------------
# pdf( "fig/res_global_abs_lim_unlim.pdf", width=6, height=5 )
  par(las=1, mar=c(4,4,1,1))
  with( adf, plot(  year, res_gpp, type='l', ylab="residual from linear trend (PgC/yr)" ) )
  with( adf, lines( year, res_gpp_lim, col='red' ) )
  legend("topright", c("GPP, not limited by soil moisture", "GPP, limited by soil moisture"), col=c("black", "red"), bty="n", lty=1 )
# dev.off()

##------------------------
## plot residuals from linear trend, dry GPP
##------------------------
  with( adf, plot(  year, -res_gpp_dry/2, type='l' ) )
  with( adf, lines( year, -res_gpp_lim_dry/2, col='red' ) )


##------------------------
## plot residuals from linear trend from normalised GPP (120)
##------------------------
  with( adf, plot(  year, res_gpp120/2, type='l' ) )
  with( adf, lines( year, res_gpp_lim120/2, col='red' ) )
