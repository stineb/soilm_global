library(ncdf4)
library(dplyr)

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

## get GPP
if (!exists("gpp")){

  # print("getting gpp ...")
  # nc <- nc_open( paste( myhome, "data/pmodel_output/Pmod_output_global_GPP_C3_s1982-01-01_e2011-12-31_r15_8.nc", sep="") )
  # gpp <- ncvar_get( nc, varid="GPP" )
  # time <- nc$dim$time$vals
  # lon <- nc$dim$longitude$vals
  # lat <- nc$dim$latitude$vals
  # units_time <- nc$dim$time$units
  # nc_close( nc )

  print("getting gpp ...")
  nc <- nc_open( "data/gpp_lim.nc" )
  gpp <- ncvar_get( nc, varid="gpp" )
  gpp_lim <- ncvar_get( nc, varid="gpp_lim" )
  time <- nc$dim$TIME$vals
  lon <- nc$dim$LONGITUDE$vals
  lat <- nc$dim$LATITUDE$vals
  units_time <- nc$dim$time$units
  nc_close( nc )

}


## get area
if (!exists("arr_area")){
  print("getting area ...")
  nc <- nc_open( "data/area.nc" )
  arr_area <- ncvar_get( nc, varid="area" )
  nc_close( nc )  
}

## get global total unlimited GPP over time
gpp_abs <- sweep( gpp, 1, arr_area, "*", check.margin=FALSE )
ggpp <- apply( gpp_abs, c(3), FUN=sum, na.rm=TRUE )

##------------------------------------------------
## construct dataframes with monthly and annual values for global GPP
##------------------------------------------------
mdf <-  data.frame( 
                  date=as.Date( time-2, origin='0001-01-01' ),
                  gpp=ggpp * 1e-15
                  ) %>%
        mutate( moy = as.numeric( format( date, format="%m" ) ),
                year = as.numeric( format( date, format="%Y" ) )
               )

adf <- mdf %>%  group_by( year ) %>% 
                summarise( gpp=sum(gpp) )


##------------------------------------------------
## get linear trends and add residuals and add to data frame
##------------------------------------------------
lm_gpp     <- lm( gpp ~ year, data=adf )
adf$iav    <- lm_gpp$residuals
adf$trend  <- lm_gpp$fitted.values - mean( lm_gpp$fitted.values )
adf$mean   <- mean( lm_gpp$fitted.values )


##------------------------------------------------
## Get mean, trend and iav for all gridcells
##------------------------------------------------
mean <- gpp[,,1]; mean[] <- NA  
trend <- gpp[,,1:30]; trend[] <- NA
iav <- gpp[,,1:30]; iav[] <- NA

mean_lim <- gpp[,,1]; mean_lim[] <- NA  
trend_lim <- gpp[,,1:30]; trend_lim[] <- NA
iav_lim <- gpp[,,1:30]; iav_lim[] <- NA

for (ilon in seq(dim(gpp)[1])){
  print(paste("ilon:", ilon))
  for (ilat in seq(dim(gpp)[2])){
    if (!is.na(gpp[ilon,ilat,1])){
      
      ## construct dataframes with monthly and annual values for global GPP
      mdf_tmp <-  data.frame( 
                        date=as.Date( time-2, origin='0001-01-01' ),
                        gpp=gpp[ilon,ilat,],
                        gpp_lim=gpp_lim[ilon,ilat,]
                        ) %>%
              mutate( moy = as.numeric( format( date, format="%m" ) ),
                      year = as.numeric( format( date, format="%Y" ) )
                     )

      adf_tmp <- mdf_tmp %>%  group_by( year ) %>% 
                      summarise( gpp=sum(gpp), gpp_lim=sum(gpp_lim) )

      if ( sum(!is.na(adf_tmp$gpp))>3 && sum(!is.na(adf_tmp$gpp_lim))>3){

        ## get linear trends and add residuals and add to data frame
        ## unlimited GPP
        lm_tmp <- lm( gpp ~ year, data=adf_tmp, na.action=na.exclude )

        adf_tmp$iav    <- residuals(lm_tmp)
        adf_tmp$trend  <- fitted(lm_tmp) - mean( fitted(lm_tmp), na.rm=TRUE )
        adf_tmp$mean   <- mean( fitted(lm_tmp), na.rm=TRUE )

        iav[ilon,ilat,] <- adf_tmp$iav
        trend[ilon,ilat,] <- adf_tmp$trend
        mean[ilon,ilat] <- adf_tmp$mean[1]

        ## limited GPP
        lm_tmp <- lm( gpp_lim ~ year, data=adf_tmp, na.action=na.exclude )

        adf_tmp$iav    <- residuals(lm_tmp)
        adf_tmp$trend  <- fitted(lm_tmp) - mean( fitted(lm_tmp), na.rm=TRUE )
        adf_tmp$mean   <- mean( fitted(lm_tmp), na.rm=TRUE )

        iav_lim[ilon,ilat,] <- adf_tmp$iav
        trend_lim[ilon,ilat,] <- adf_tmp$trend
        mean_lim[ilon,ilat] <- adf_tmp$mean[1]        

      }

    }
  }
}

## Get IAV variance fields
iav_var_rel <- apply( iav, c(1,2), FUN=var ) / mean
iav_var_lim_rel <- apply( iav_lim, c(1,2), FUN=var ) / mean_lim
iav_change <- iav_var_lim_rel / iav_var_rel 

##------------------------
## Plot map for IAV increase
##------------------------
  magn <- 4
  ncols <- 2
  nrows <- 1
  widths <- rep(1.6*magn,ncols)
  widths[2] <- 0.25*widths[1]
  heights <- rep(magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  ylim <- c(-60,85)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)

  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )

  pdf( "fig/map_iav_increase.pdf", width=sum(widths), height=sum(heights) )

    panel <- layout(
              order,
              widths=widths,
              heights=heights,
              TRUE
              )
    # layout.show( panel )

    ## Color key
    # par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    color <- c( "royalblue4", "royalblue1", "wheat", "tomato", "tomato4" )
    # lev <- c( 0, 0.5, 0.7, 0.8, 0.9, 0.95, 1, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5, 2, 999 )
    # lev <- seq(-2,2,0.2)
    lev <- c( 0, 0.2, 0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.5, 2, 3, 999 )
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=FALSE, maxval=1000 )

    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    image(
            seq(-179.75, 179.75, 0.5), seq(-89.75, 89.75, 0.5), 
            iav_change,
            ylim=c(-60,85), 
            # zlim=range(lev), 
            yaxt="n", xaxt="n",
            col=out.mycolorbar$colors, breaks=out.mycolorbar$margins,
            xlab="", ylab=""
            )
    map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )

    axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
    axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )

    axis( 4, at=lat.labels, lab=F, lwd=1.5 )
    axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )

    axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
    axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )

    axis( 3, at=lon.labels, lab=F, lwd=1.5 )
    axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )

    ## Color key
    par( mar=c(3,3,1,1),xaxs="i", yaxs="i",las=1)
    out.mycolorbar <- mycolorbar( color, lev, orient="v", plot=TRUE, maxval=1 )

  dev.off()  

