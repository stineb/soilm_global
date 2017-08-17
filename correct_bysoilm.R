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

# load( "/alphadata01/bstocker/sofun/utils_sofun/analysis_sofun/fluxnet2015/data/nice_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'

# load( "/alphadata01/bstocker/sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

# ## add date and MOY to dataframe nice_agg
# nice_agg <- nice_agg %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
# nice_agg <- nice_agg %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

# ## add alpha to nice_agg
# nice_agg <- nice_agg %>% mutate( alpha = aet_pmodel / pet_pmodel )

# ## aggregate nice_agg to monthly values
# mdf <- nice_agg %>% group_by( mysitename, year, moy ) %>% summarise( flue = mean( flue, na.rm=TRUE ), soilm = mean( soilm_mean, na.rm=TRUE ) )

# ## Bin values and get mean fLUE for soil moisture < 0.25 for each site (:= flue0)
# intervals <- seq(0, 1, 0.25) #c( 0, 0.15, 0.85, 1.0 )
# mid <- intervals[1:(length(intervals)-1)] + (intervals[2]-intervals[1])/2
# mdf$ininterval <- NULL
# mdf <- mdf %>% mutate( ininterval = cut( soilm , breaks = intervals ) ) %>% group_by( mysitename, ininterval )
# mdf_agg <- mdf %>% dplyr::summarise( flue0=mean( flue, na.rm=TRUE ) ) %>% complete( ininterval, fill = list( flue0 = NA ) )# %>% dplyr::select( median )
# mdf_flue0 <- mdf_agg %>% dplyr::filter( ininterval=="(0,0.25]" )
# # flue_vs_soilm[ jdx, ] <- unlist(tmp)[1:nintervals]

# ## Merge alpha values into this dataframe
# mdf_flue0 <- mdf_flue0 %>% left_join( df_alpha, by="mysitename" )

# ## get beta factor for quadratic fit 
# ## flue0 = beta * (x0 - x1)^2 + 1
# x0 <- mid[1]
# x1 <- 0.75
# mdf_flue0 <- mdf_flue0 %>% mutate( beta = (flue0 - 1) / (x0 - x1)^2 )
 
# ##------------------------------------
# ## plot flue0 vs. alpha
# ##------------------------------------
#   with( mdf_flue0, plot( alpha, flue0, pch=16 ) )
#   linmod <- lm( flue0 ~ alpha, data=mdf_flue0 )
#   abline( linmod, col="black" )

#   text( 0.3, 1.2, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
#   cf <- coef(linmod) %>% round( 2 )
#   eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
#   text( 0.3, 1.1, eq, adj=0.0 )


# ##------------------------------------
# ## for each site, plot fit
# ##------------------------------------
#   stress_quad <- function( x, beta ){
#     x0 <- 0.125
#     x1 <- 0.75
#     if (x>x1) {
#       y <- 1
#     } else {
#       y <- beta * ( x - x1 )^2 + 1.0
#     }
#     return( y )
#   }

#   mycurve <- function( func, from, to, col='black', add=FALSE, lwd=1, lty=1 ){
#     range_x <- seq( from, to, by=(to-from)/100 )
#     range_y <- sapply( range_x, func )
#     if (add){
#       lines( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
#     } else {
#       plot( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
#     }
#   }

#   for ( sitename in unique(mdf$mysitename) ){
#     with( dplyr::filter(mdf, mysitename==sitename),
#           plot(
#                 soilm, 
#                 flue,
#                 pch=16,
#                 col=rgb(0,0,1,0.6),
#                 main=sitename,
#                 ylim=c(0,1.4),
#                 xlim=c(0,1),
#                 xlab="soil water content (fraction)",
#                 ylab="fLUE"            
#                 )
#           )
#     mybeta <- dplyr::filter( mdf_flue0, mysitename==sitename)$beta
#     if (!is.na(mybeta)) mycurve( function(x) stress_quad( x, mybeta ), from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
#   }

# ##------------------------------------
# ## plot scatterplot of flue vs soilm 
# ##------------------------------------
# # par( las=1, mar=c(4,4,1,1), xaxs="i", yaxs="i" )
#   par( las=1 )
#   with( 
#         mdf,
#         heatscatter( 
#                     soilm, 
#                     flue,
#                     ylim=c(0,1.4),
#                     xlim=c(0,1),
#                     main="",
#                     xlab="soil water content (fraction)",
#                     ylab="fLUE", 
#                     cex.lab=1.2
#                     ) 

#       )
#   box( lwd=1.5 )
#   axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
#   axis( 2, lwd=1.5, cex.lab=2 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
#   axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
#   axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
#   abline( h=1.0, lwd=0.5 )

#   for ( sitename in unique(mdf$mysitename) ){
#     mybeta <- dplyr::filter( mdf_flue0, mysitename==sitename)$beta
#     if (!is.na(mybeta)) mycurve( function(x) stress_quad( x, mybeta ), from=0.0, to=1.0, col='black', add=TRUE, lwd=1 )
#   }

#   legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )
