library(ncdf4)
library(dplyr)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

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

flue <- soilm
flue[] <- NA
for (ilon in seq(dim(aalpha)[1])){
  for (ilat in seq(dim(aalpha)[2])){
    if (!is.na(aalpha[ilon,ilat])){
      flue[ilon,ilat,] <- sapply( soilm[ilon,ilat,], FUN = function(x) calc_flue( x, beta[ilon,ilat] ) )
    }
  }
}

# ## write to file
# cdf.write( flue, "flue", 
#            lon, lat,
#            time = time,
#            make.tdim = TRUE,
#            units_time = units_time,
#            long_name_var1 = "fractional reduction of LUE due to soil moisture",
#            units_var1 = "fraction",
#            filnam = "data/flue.nc" 
#            )


## get global total unlimited GPP over time
gpp_abs <- sweep( gpp, 1, arr_area, "*", check.margin=FALSE )
ggpp <- apply( gpp_abs, c(3), FUN=sum, na.rm=TRUE )

## get global total unlimited GPP over time
gpp_lim <- gpp * flue
gpp_lim_abs <- sweep( gpp_lim, 1, arr_area, "*", check.margin=FALSE )
ggpp_lim <- apply( gpp_lim_abs, c(3), FUN=sum, na.rm=TRUE )

## write to file
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

## construct dataframes with monthly and annual values
mdf <-  data.frame( 
                  date=as.Date( time-2, origin='0001-01-01' ),
                  gpp=ggpp * 1e-15,
                  gpp_lim=ggpp_lim * 1e-15
                  ) %>%
        mutate( moy = as.numeric( format( date, format="%m" ) ),
                year = as.numeric( format( date, format="%Y" ) )
               )

adf <- mdf %>% group_by( year ) %>% summarise( gpp=sum(gpp), gpp_lim=sum(gpp_lim) ) %>% mutate( gpp_inc=gpp/gpp[1], gpp_inc_lim=gpp_lim/gpp_lim[1] )

## get linear trends and add residuals to data frame
lm_gpp     <- lm( gpp ~ year, data=adf )
lm_gpp_lim <- lm( gpp_lim ~ year, data=adf )
adf$res_gpp <- lm_gpp$residuals
adf$res_gpp_lim <- lm_gpp_lim$residuals

lm_gpp_inc     <- lm( gpp_inc ~ year, data=adf )
lm_gpp_inc_lim <- lm( gpp_inc_lim ~ year, data=adf )
adf$res_gpp_inc <- lm_gpp_inc$residuals
adf$res_gpp_inc_lim <- lm_gpp_inc_lim$residuals

## plot time series, absolute
  par(las=1)
  with( adf, plot(year, gpp, type="l", ylim=c(130,180), ylab="GPP (PgC/yr)" ) )
  abline( lm_gpp, col="black", lty=2 )
  cf <- coef(lm_gpp) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-1, adj=0.1 )

  with( adf, lines(year, gpp_lim, col="red" ) )
  lm_gpp_lim <- lm( gpp_lim ~ year, data=adf )
  abline( lm_gpp_lim, col="red", lty=2 )
  cf <- coef(lm_gpp_lim) %>% round( 2 )
  eq <- paste0( "slope = ", abs(cf[2]), " PgC/yr " )
  mtext( eq, line=-2, adj=0.1, col='red' )

## plot time series, relative increase
  par(las=1)
  with( adf, plot(year, gpp_inc, type="l", ylim=c(0.95,1.2), ylab="change in GPP (frac., rel. to 1982)" ) )
  abline( lm_gpp_inc, col="black", lty=2 )
  cf <- coef(lm_gpp_inc) %>% round( 4 )
  eq <- paste0( "slope = ", abs(cf[2]), " yr-1 " )
  mtext( eq, line=-1, adj=0.1 )

  with( adf, lines(year, gpp_inc_lim, col="red" ) )
  lm_gpp_inc_lim <- lm( gpp_inc_lim ~ year, data=adf )
  abline( lm_gpp_inc_lim, col="red", lty=2 )
  cf <- coef(lm_gpp_inc_lim) %>% round( 4 )
  eq <- paste0( "slope = ", abs(cf[2]), " yr-1 " )
  mtext( eq, line=-2, adj=0.1, col='red' )

## plot flue effect
with( adf, plot( year, gpp_lim / gpp, type="l" ) )

## plot residuals from linear trend
with( adf, plot( year, -res_gpp*1e-15/2, type='l' ) )
with( adf, lines( year, -res_gpp_lim*1e-15/2, col='red' ) )
