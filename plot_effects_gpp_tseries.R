library(ncdf4)
library(RColorBrewer)
library(dplyr)
source("plot_map.R")
source("../utilities/add_alpha.R")


overwrite <- FALSE

filn <- "data/gpp_glob_tseries.Rdata"

if (!file.exists(filn)||overwrite){
  ##------------------------------------------------------------------------
  ## GPP loss
  ##------------------------------------------------------------------------
  fil_s0  <- "gpp_pmodel_s0_GLOB.nc"
  fil_s1a <- "gpp_pmodel_s1a_GLOB.nc"
  fil_s1b <- "gpp_pmodel_s1b_GLOB.nc"
  fil_s1c <- "gpp_pmodel_s1c_GLOB.nc"

  dir <- paste0( myhome, "/data/pmodel_fortran_output/v2/")

  nc <- nc_open( paste0( dir, fil_s0 ) )
  gpp_s0 <- ncvar_get( nc, varid="gpp" )
  time <- nc$dim$time$vals
  nc_close(nc)

  ## s1a
  nc <- nc_open( paste0( dir, fil_s1a ) )
  gpp_s1a <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)

  ## s1b
  nc <- nc_open( paste0( dir, fil_s1b ) )
  gpp_s1b <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)

  ## s1c
  nc <- nc_open( paste0( dir, fil_s1c ) )
  gpp_s1c <- ncvar_get( nc, varid="gpp" )
  nc_close(nc)    

  ## modis
  nc <- nc_open( paste0( myhome, "/data/gpp_modis/gpp_modis_GLOB.nc" ) )
  gpp_modis <- ncvar_get( nc, varid="gpp" )
  year_modis <- 2000:2015
  df_modis <- tibble( year=year_modis, gpp_modis=gpp_modis)
  nc_close(nc)    

  ## vpm
  nc <- nc_open( paste0( myhome, "/data/gpp_vpm/gpp_vpm_GLOB.nc" ) )
  gpp_vpm <- ncvar_get( nc, varid="gpp" )
  year_vpm <- 2000:2016
  df_vpm <- tibble( year=year_vpm, gpp_vpm=gpp_vpm)
  nc_close(nc)    

  ## bess
  nc <- nc_open( paste0( myhome, "/data/gpp_bess/gpp_bess_GLOB.nc" ) )
  gpp_bess <- ncvar_get( nc, varid="gpp" )
  year_bess <- 2001:2015
  df_bess <- tibble( year=year_bess, gpp_bess=gpp_bess)
  nc_close(nc)   

  ## mte
  nc <- nc_open( paste0( myhome, "/data/gpp_mte/gpp_mte_GLOB.nc" ) )
  gpp_mte <- ncvar_get( nc, varid="gpp" )
  year_mte <- 1982:2011
  df_mte <- tibble( year=year_mte, gpp_mte=gpp_mte)
  nc_close(nc)        

  year <- 1982:2016
  
  ## create data frame
  df <- tibble( year=year, gpp_s0=gpp_s0, gpp_s1a=gpp_s1a, gpp_s1b=gpp_s1b, gpp_s1c=gpp_s1c ) %>%
        left_join( df_modis, by = "year") %>%
        left_join( df_vpm, by = "year") %>%
        left_join( df_bess, by = "year") %>%
        left_join( df_mte, by = "year")
        
  save( df, file = filn)

} else {

  load( filn )

}

##-----------------------------------------------------
## Get additional variables
##-----------------------------------------------------
## difference due to soil moisture effects
df <- df %>% mutate( diff_s1a = gpp_s1a - gpp_s0, diff_s1b = gpp_s1b - gpp_s0, diff_s1c = gpp_s1c - gpp_s0 )

## linear trends
linmod_s0  <- lm( gpp_s0  ~ year, data=df )
linmod_s1a <- lm( gpp_s1a ~ year, data=df )
linmod_s1b <- lm( gpp_s1b ~ year, data=df )
linmod_s1c <- lm( gpp_s1c ~ year, data=df )
  
linmod_s0_rel  <- lm( 100*(gpp_s0 / mean(gpp_s0) - 1)  ~ year, data=df )
linmod_s1a_rel <- lm( 100*(gpp_s1a / mean(gpp_s1a) - 1) ~ year, data=df )
linmod_s1b_rel <- lm( 100*(gpp_s1b / mean(gpp_s1b) - 1) ~ year, data=df )
linmod_s1c_rel <- lm( 100*(gpp_s1c / mean(gpp_s1c) - 1) ~ year, data=df )

linmod_s1a_releff <- lm( -100*(1-gpp_s1a/gpp_s0) ~ year, data=df )
linmod_s1b_releff <- lm( -100*(1-gpp_s1b/gpp_s0) ~ year, data=df )
linmod_s1c_releff <- lm( -100*(1-gpp_s1c/gpp_s0) ~ year, data=df )

##-----------------------------------------------------
## GPP time series
##-----------------------------------------------------
plot_gpp_global_tseries <- function( df, filn=NA ){
  if (!is.na(filn)) pdf(, width = 7, height = 6 )
    par(las=1)
    with( df, plot( year, gpp_s1b, type="l", ylim=c(100,160), lty=2, xlab="year", ylab=expression( paste("Global GPP (PgC yr"^-1, ")" ) ) ) )
    with( df,  polygon( c( year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), border = NA, col=rgb(0,0,0,0.3) ) )

    with( df, lines( year, gpp_s0 ) )
    with( df, lines( year, gpp_modis, col="red" ) )
    with( df, lines( year, gpp_vpm, col="magenta" ) )
    with( df, lines( year, gpp_bess, col="blue" ) )
    with( df, lines( year, gpp_mte, col="green" ) )
    axis(4, labels = FALSE)

    legend("topleft", c( "P-model", "P-model, corrected (IV)", "MODIS", "VPM", "BESS", "MTE" ), bty = "n", lty = c(1,2,1,1,1,1), col=c("black", "black", "red", "magenta", "blue", "green") )
  if (!is.na(filn)) dev.off()  
}
plot_gpp_global_tseries( df, filn="fig/gpp_global_tseries.pdf" )

##-----------------------------------------------------
## Absolute trends with P-model output
##-----------------------------------------------------
plot_gpp_global_tseries_trend_pmodel <- function( df, filn=NA ){
  if (!is.na(filn)) pdf(filn, width = 7, height = 6 )
    par(las=1, xaxs="i", yaxs="i")
    with( df, plot( year, gpp_s1b, type="l", ylim=c(100,160), col="red", lty=1, xlab="year", ylab=expression( paste("Global GPP (PgC yr"^-1, ")" ) ) ) )
    abline( linmod_s1b, col="red" )
    cf <- coef(linmod_s1b) %>% round( 2 )
    eq <- paste0( "s1b: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-2, adj=0.05, col="red" )

    with( df, polygon( c(year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), col=add_alpha("red",0.5), border=NA ) )
    abline( linmod_s1a, col=add_alpha("red",0.5) )

    abline( linmod_s1c, col=add_alpha("red",0.5) )

    with( df, lines( year, gpp_s0 ) )
    abline( linmod_s0 )
    cf <- coef(linmod_s0) %>% round( 2 )
    eq <- paste0( "s0: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-1, adj=0.05, col="black" )
    axis(4, labels = FALSE)
    legend( "bottomleft", c("s0", "s1b"), col=c("black", "red"), bty = "n", lty=1, lwd=1 )
  if (!is.na(filn)) dev.off()  
} 
plot_gpp_global_tseries_trend_pmodel( df, filn="fig/gpp_global_tseries_trend_pmodel.pdf" )

##-----------------------------------------------------
## Relative trends with P-model output
##-----------------------------------------------------
plot_gpp_global_tseries_reltrend_pmodel <- function( df, filn=NA ){
  if (!is.na(filn)) pdf(filn, width = 7, height = 6 )
    par(las=1, xaxs="i", yaxs="i")
    with( df, plot( year, 100*(gpp_s1b/mean(gpp_s1b)-1), type="l", col="red", lty=1, xlab="year", ylab=expression( paste("Global GPP (% change)" ) ), ylim=c(-10,8) ) )
    abline( linmod_s1b_rel, col="red" )
    cf <- coef(linmod_s1b_rel) %>% round( 4 )
    eq <- paste0( "s0: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-3, adj=0.05, col="red" )

    with( df, polygon( c(year, rev(year)), c(100*(gpp_s1a/mean(gpp_s1a)-1), rev(100*(gpp_s1c/mean(gpp_s1c)-1))), col=add_alpha("red", 0.5), border = NA ) )
    abline( linmod_s1a_rel, col=add_alpha("red", 0.5) )
    cf <- coef(linmod_s1a_rel) %>% round( 4 )
    eq <- paste0( "s1a: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-2, adj=0.05, col=add_alpha("red", 0.5) )

    abline( linmod_s1c_rel, col=add_alpha("red", 0.5) )
    cf <- coef(linmod_s1c_rel) %>% round( 4 )
    eq <- paste0( "s1b: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-4, adj=0.05, col=add_alpha("red", 0.5) )

    with( df, lines( year, 100*(gpp_s0/mean(gpp_s0)-1) ) )
    abline( linmod_s0_rel )
    cf <- coef(linmod_s0_rel) %>% round( 4 )
    eq <- paste0( "s1c: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-1, adj=0.05, col="black" )
    axis(4, labels = FALSE)
    legend( "bottomright", c("s0", "s1b"), col=c("black", "red"), bty = "n", lty=1, lwd=1 )
  if (!is.na(filn)) dev.off()
}
plot_gpp_global_tseries_reltrend_pmodel( df, filn="fig/gpp_global_tseries_reltrend_pmodel.pdf" )

##-----------------------------------------------------
## Trends of soil moisture effect P-model output
##-----------------------------------------------------
plot_gpp_global_tseries_relefftrend_pmodel <- function( df, filn=NA ){
  if (!is.na(filn)) pdf(filn, width = 7, height = 6 )
    par(las=1, xaxs="i", yaxs="i")
    with( df, plot( year, -100*(1-gpp_s1b/gpp_s0), type="l", col="red", lty=1, xlab="year", ylim=c(-22, -8), ylab=expression( paste("reduction in global GPP (%)" ) ) ) )
    abline( linmod_s1b_releff, col="red" )
    cf <- coef(linmod_s1b_releff) %>% round( 4 )
    eq <- paste0( "s1b: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-3, adj=0.05, col="red" )

    with( df, polygon( c(year, rev(year)), c( -100*(1-gpp_s1c/gpp_s0), rev( -100*(1-gpp_s1a/gpp_s0)) ), col=add_alpha("red", 0.5), border = NA ) )
    abline( linmod_s1a_releff, col=add_alpha("red", 0.5) )
    cf <- coef(linmod_s1a_releff) %>% round( 4 )
    eq <- paste0( "s1a: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-2, adj=0.05, col=add_alpha("red", 0.5) )

    abline( linmod_s1c_releff, col=add_alpha("red", 0.5) )
    cf <- coef(linmod_s1c_releff) %>% round( 4 )
    eq <- paste0( "s1c: slope = ", ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]) )
    mtext( eq, line=-4, adj=0.05, col=add_alpha("red", 0.5) )

    axis(4, labels = FALSE)

  if (!is.na(filn)) dev.off()
}
plot_gpp_global_tseries_relefftrend_pmodel( df, filn="fig/gpp_global_tseries_relefftrend_pmodel.pdf" )


get_numbers_effects_gpp_tseries <- function( df ){
  with( df, print( paste0( "Percent decrease in mean GPP across years 1982-2016: ", 
    format( 100*(1 - mean(gpp_s1b)/mean(gpp_s0)), digits=4), "% (", 
    format( 100*(1 - mean(gpp_s1a)/mean(gpp_s0)), digits=4), " - ", 
    format( 100*(1 - mean(gpp_s1c)/mean(gpp_s0)), digits=4), "%)" ) ) )
} 


