library(dplyr)
library(tidyr)
library(LSD)
library(minpack.lm)

load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'

load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

## get flue0 : maximum fLUE reduction at low soil moisture
calc_flue0 <- function( alpha ){ 
  ## equation from correct_bysoilm.R
  y <- 0.1 + 0.92 * alpha
  y <- min( y, 1.0 )
  y <- max( y, 0.0 )
  return( y )
}

calc_beta <- function( flue0 ){
  x1 <- 0.75
  x0 <- 0.125
  beta <- (flue0 - 1) / (x0 - x1)^2
}

calc_flue_est <- function( x, flue0 ){
  if (is.na(flue0)){
    y <- 1.0
  } else if ( flue0>1.0 ){
    y <- 1
  } else {
    beta <- calc_beta( flue0 )
    x1 <- 0.75
    if (x>x1) {
      y <- 1
    } else {
      y <- beta * ( x - x1 )^2 + 1.0
    }    
  }
  return( y )
}

calc_flue_est_alpha <- function( x, alpha, apar, bpar, cpar, dpar ){
  # alpha <- min( 1.0, alpha )
  flue0 <- apar + bpar * alpha
  # flue0 <- min( 1.0, max( 0.0, flue0 ) )
  beta  <- (flue0 - 1.0) / (cpar - dpar)^2
  flue_est <- beta * (x - dpar)^2 + 1
  return( flue_est )
}

mycurve <- function( func, from, to, col='black', add=FALSE, lwd=1, lty=1 ){
  range_x <- seq( from, to, by=(to-from)/100 )
  range_y <- sapply( range_x, func )
  if (add){
    lines( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
  } else {
    plot( range_x, range_y, type="l", col=col, lwd=lwd, lty=lty )
  }
}

## add date and MOY to dataframe nice_agg
ddf <- nice_agg %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
ddf <- ddf      %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

## aggregate nice_agg to monthly values
mdf <- ddf %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_mean = mean( soilm_mean, na.rm=TRUE ) )

##------------------------------------
## Estimate fLUE using minimum fLUE ~ alpha relationship determined independently
##------------------------------------
## Bin values and get mean fLUE for soil moisture < 0.25 for each site (:= flue0)
## monthly data
intervals <- seq(0, 1, 0.25) #c( 0, 0.15, 0.85, 1.0 )
mid <- intervals[1:(length(intervals)-1)] + (intervals[2]-intervals[1])/2
mdf$ininterval <- NULL
mdf <- mdf %>% mutate( ininterval = cut( soilm_mean , breaks = intervals ) ) %>% group_by( mysitename, ininterval )
mdf_agg <- mdf %>% dplyr::summarise( flue0=mean( fvar, na.rm=TRUE ) ) %>% complete( ininterval, fill = list( flue0 = NA ) )# %>% dplyr::select( median )
mdf_flue0 <- mdf_agg %>% dplyr::filter( ininterval=="(0,0.25]" )

## Merge alpha values into this dataframe
mdf_flue0 <- mdf_flue0 %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

## get beta factor for quadratic fit 
## flue0 = beta * (x0 - x1)^2 + 1
mdf_flue0 <- mdf_flue0 %>% mutate( flue0_est = calc_flue0( meanalpha ) )
 
## Estimate fLUE based on soil moisture and alpha
mdf <- mdf %>%  left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" ) %>%
                left_join( select( mdf_flue0, flue0, flue0_est ), by="mysitename" ) %>%
                mutate( flue_est = calc_flue_est( soilm_mean, flue0_est ) )

## Daily data:
ddf$ininterval <- NULL
ddf <- ddf %>% mutate( ininterval = cut( soilm_mean , breaks = intervals ) ) %>% group_by( mysitename, ininterval )
ddf_agg <- ddf %>% dplyr::summarise( flue0=mean( fvar, na.rm=TRUE ) ) %>% complete( ininterval, fill = list( flue0 = NA ) )# %>% dplyr::select( median )
ddf_flue0 <- ddf_agg %>% dplyr::filter( ininterval=="(0,0.25]" )

## Merge alpha values into this dataframe
ddf_flue0 <- ddf_flue0 %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

## get beta factor for quadratic fit 
## flue0 = beta * (x0 - x1)^2 + 1
ddf_flue0 <- ddf_flue0 %>% mutate( flue0_est = calc_flue0( meanalpha ) )
 
## Estimate fLUE based on soil moisture and alpha
ddf <- ddf %>%  left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" ) %>%
                left_join( select( mdf_flue0, flue0, flue0_est ), by="mysitename" ) %>%
                mutate( flue_est = calc_flue_est_alpha( soilm_mean, meanalpha, 0.11, 0.89, 0.125, 0.75 ) )

##------------------------------------
## Estimate fLUE using non-linear least squares on quadratic function
##------------------------------------
## Fit by medians in bis
# df_tmp <- data.frame( xvals=xvals, yvals=yvals )
# df_tmp <- filter( ddf, soilm_mean, meanalpha, fvar )
stressfit <- try( 
                  nlsLM( 
                        fvar ~ calc_flue_est_alpha( soilm_mean, meanalpha, apar, bpar, cpar, dpar ),
                        data=ddf,
                        start=list( apar=0.11, bpar=0.89, cpar=0.125, dpar=0.75 ),
                        lower=c( apar=0.001, bpar=0.5, cpar=0.001, dpar=0.3 ),
                        upper=c( apar=0.3, bpar=1.5, cpar=0.3, dpar=1.0 ),
                        algorithm="port"
                        ) 
                  )

print("fitted coefficients are:")
print(coef(stressfit))

ddf <- ddf %>% mutate( flue_est_nls = calc_flue_est_alpha( soilm_mean, meanalpha, coef(stressfit)[[ "apar" ]], coef(stressfit)[[ "bpar" ]], coef(stressfit)[[ "cpar" ]], coef(stressfit)[[ "dpar" ]] ) )

##------------------------------------
## plot estimated and "actual" fLUE
##------------------------------------
  with( ddf, plot( fvar, flue_est, pch=16, col=rgb(0,0,0,0.1)) )
  lines( c(-100,100), c(-100,100), lty=3, col="red" )

  with( ddf, plot( fvar, flue_est_nls, pch=16, col=rgb(0,0,0,0.1)) )
  lines( c(-100,100), c(-100,100), lty=3, col="red" )  

  with( filter( ddf, mysitename=="US-Ton"), plot( year_dec, fvar, type="l" ) )
  with( filter( ddf, mysitename=="US-Ton"), lines( year_dec, flue_est, col="red" ) )
  with( filter( ddf, mysitename=="US-Ton"), lines( year_dec, flue_est_nls, col="green" ) )

##------------------------------------
## plot flue0 vs. alpha
##------------------------------------
## Monthly data:
  par( las=1 )
  with( mdf_flue0, plot( meanalpha, flue0, pch=16, xlab="AET/PET", ylab=expression(paste("fLUE"[0])) ) )
  linmod <- lm( flue0 ~ meanalpha, data=mdf_flue0 )
  abline( linmod, col="black" )
  abline( a=coef(stressfit)[[ "apar" ]], b=coef(stressfit)[[ "bpar" ]], col="red" )

  text( 0.3, 1.2, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
  cf <- coef(linmod) %>% round( 2 )
  eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  text( 0.3, 1.1, eq, adj=0.0 )
  eqfit <- paste0( "y = ", format( coef(stressfit)[[ "apar" ]], digits=2 ), ifelse(sign(coef(stressfit)[[ "bpar" ]])==1, " + ", " - "), format( abs(coef(stressfit)[[ "bpar" ]]), digits=2 ), " x " )
  text( 0.3, 1.0, eqfit, adj=0.0, col="red" )
  
  

  ## some sites have fLUE0 = NA - these are mostly the ones with high AET/PET
  boxplot( alpha ~ is.na(flue0), data=mdf_flue0 )
  
## Daily data:
  par( las=1 )
  with( ddf_flue0, plot( meanalpha, flue0, pch=16, xlab="AET/PET", ylab=expression(paste("fLUE"[0])) ) )
  linmod <- lm( flue0 ~ meanalpha, data=ddf_flue0 )
  abline( linmod, col="black" )
  abline( a=coef(stressfit)[[ "apar" ]], b=coef(stressfit)[[ "bpar" ]], col="red" )

  text( 0.3, 1.2, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
  cf <- coef(linmod) %>% round( 2 )
  eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  text( 0.3, 1.1, eq, adj=0.0 )
  eqfit <- paste0( "y = ", format( coef(stressfit)[[ "apar" ]], digits=2 ), ifelse(sign(coef(stressfit)[[ "bpar" ]])==1, " + ", " - "), format( abs(coef(stressfit)[[ "bpar" ]]), digits=2 ), " x " )
  text( 0.3, 1.0, eqfit, adj=0.0, col="red" )

  ## some sites have fLUE0 = NA - these are mostly the ones with high AET/PET
  boxplot( alpha ~ is.na(flue0), data=ddf_flue0 )

##------------------------------------
## for each site, plot fit
##------------------------------------
  ## Monthly data
  ## Plot relationship between soil moisture and fLUE for each site individually
  for ( sitename in unique(mdf$mysitename) ){
    with( dplyr::filter(mdf, mysitename==sitename),
          plot(
                soilm_mean, 
                fvar,
                pch=16,
                col=rgb(0,0,1,0.6),
                main=sitename,
                ylim=c(0,1.4),
                xlim=c(0,1),
                xlab="soil water content (fraction)",
                ylab="fLUE"            
                )
          )
    abline( h=1.0, lty=3 )
    mycurve( function(x) calc_flue_est( x, dplyr::filter( mdf_flue0, mysitename==sitename )$flue0 ),     from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
    mycurve( function(x) calc_flue_est( x, dplyr::filter( mdf_flue0, mysitename==sitename )$flue0_est ), from=0.0, to=1.0, col='orange', add=TRUE, lwd=2, lty=2 )

    mycurve( function(x) calc_flue_est_alpha( x, meanalpha, coef(stressfit)[[ "apar" ]], coef(stressfit)[[ "bpar" ]], coef(stressfit)[[ "cpar" ]], coef(stressfit)[[ "dpar" ]] ), from=0.0, to=1.0, col='green', add=TRUE, lwd=2 )

  }

  ## Daily data
  ## Plot relationship between soil moisture and fLUE for each site individually
  for ( sitename in unique(ddf$mysitename) ){
    with( dplyr::filter(ddf, mysitename==sitename),
          plot(
                soilm_mean, 
                fvar,
                pch=16,
                col=rgb(0,0,1,0.6),
                main=sitename,
                ylim=c(0,1.4),
                xlim=c(0,1),
                xlab="soil water content (fraction)",
                ylab="fLUE"            
                )
          )
    abline( h=1.0, lty=3 )
    mycurve( function(x) calc_flue_est( x, dplyr::filter( ddf_flue0, mysitename==sitename )$flue0 ),     from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
    mycurve( function(x) calc_flue_est( x, dplyr::filter( ddf_flue0, mysitename==sitename )$flue0_est ), from=0.0, to=1.0, col='orange', add=TRUE, lwd=2, lty=2 )
  }

# ##------------------------------------
# ## plot scatterplot of flue vs soilm_mean 
# ##------------------------------------
# # par( las=1, mar=c(4,4,1,1), xaxs="i", yaxs="i" )
#   par( las=1 )
#   with( 
#         mdf,
#         heatscatter( 
#                     soilm_mean, 
#                     flue_est,
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
#     if (!is.na(mybeta)) mycurve( function(x) calc_flue_est( x, mybeta ), from=0.0, to=1.0, col='black', add=TRUE, lwd=1 )
#   }

#   legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )
