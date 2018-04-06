library(ncdf4)
library(RColorBrewer)
library(dplyr)
source("plot_map.R")
source("~/.Rprofile")

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
  
  ##-----------------------------------------------------
  ## GPP time series
  ##-----------------------------------------------------
  pdf("fig/gpp_global_tseries.pdf", width = 7, height = 6 )
    par(las=1)
    with( df, plot( year, gpp_s1b, type="l", ylim=c(100,160), lty=2, xlab="year", ylab=expression( paste("global GPP (PgC yr"^-1, ")" ) ) ) )
    with( df,  polygon( c( year, rev(year)), c(gpp_s1a, rev(gpp_s1c)), border = NA, col=rgb(0,0,0,0.3) ) )

    with( df, lines( year, gpp_s0 ) )
    with( df, lines( year, gpp_modis, col="red" ) )
    with( df, lines( year, gpp_vpm, col="magenta" ) )
    with( df, lines( year, gpp_bess, col="blue" ) )
    with( df, lines( year, gpp_mte, col="green" ) )

    legend("topleft", c( "P-model", "P-model, corrected (IV)", "MODIS", "VPM", "BESS", "MTE" ), bty = "n", lty = c(1,2,1,1,1,1), col=c("black", "black", "red", "magenta", "blue", "green") )
  dev.off()

