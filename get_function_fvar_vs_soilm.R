load( "/alphadata01/bstocker/sofun/utils_sofun/analysis_sofun/fluxnet2015/data/nice_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'

load( "/alphadata01/bstocker/sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'

## add date and MOY to dataframe nice_agg
nice_agg <- nice_agg %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
nice_agg <- nice_agg %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

## add alpha to nice_agg
nice_agg <- nice_agg %>% mutate( alpha = aet_pmodel / pet_pmodel )

## aggregate nice_agg to monthly values
mdf <- nice_agg %>% group_by( mysitename, year, moy ) %>% summarise( flue = mean( flue, na.rm=TRUE ), soilm = mean( soilm_mean, na.rm=TRUE ) )

## Bin values and get mean fLUE for soil moisture < 0.25 for each site (:= flue0)
intervals <- seq(0, 1, 0.25) #c( 0, 0.15, 0.85, 1.0 )
mid <- intervals[1:(length(intervals)-1)] + (intervals[2]-intervals[1])/2
mdf$ininterval <- NULL
mdf <- mdf %>% mutate( ininterval = cut( soilm , breaks = intervals ) ) %>% group_by( mysitename, ininterval )
mdf_agg <- mdf %>% dplyr::summarise( flue0=mean( flue, na.rm=TRUE ) ) %>% complete( ininterval, fill = list( flue0 = NA ) )# %>% dplyr::select( median )
mdf_flue0 <- mdf_agg %>% dplyr::filter( ininterval=="(0,0.25]" )
# flue_vs_soilm[ jdx, ] <- unlist(tmp)[1:nintervals]

## Merge alpha values into this dataframe
mdf_flue0 <- mdf_flue0 %>% left_join( df_alpha, by="mysitename" )

## get beta factor for quadratic fit 
## flue0 = beta * (x0 - x1)^2 + 1
x0 <- mid[1]
x1 <- 0.75
mdf_flue0 <- mdf_flue0 %>% mutate( beta = (flue0 - 1) / (x0 - x1)^2 )
 
##------------------------------------
## plot flue0 vs. alpha
##------------------------------------
  with( mdf_flue0, plot( alpha, flue0, pch=16 ) )
  linmod <- lm( flue0 ~ alpha, data=mdf_flue0 )
  abline( linmod, col="black" )

  text( 0.3, 1.2, bquote( italic(R)^2 == .(format( summary( linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
  cf <- coef(linmod) %>% round( 2 )
  eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  text( 0.3, 1.1, eq, adj=0.0 )


##------------------------------------
## for each site, plot fit
##------------------------------------
  stress_quad <- function( x, beta ){
    x1 <- 0.75
    if (x>x1) {
      y <- 1
    } else {
      y <- beta * ( x - x1 )^2 + 1.0
    }
    return( y )
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

  for ( sitename in unique(mdf$mysitename) ){
    with( dplyr::filter(mdf, mysitename==sitename),
          plot(
                soilm, 
                flue,
                pch=16,
                col=rgb(0,0,1,0.6),
                main=sitename,
                ylim=c(0,1.4),
                xlim=c(0,1),
                xlab="soil water content (fraction)",
                ylab="fLUE"            
                )
          )
    mybeta <- dplyr::filter( mdf_flue0, mysitename==sitename)$beta
    if (!is.na(mybeta)) mycurve( function(x) stress_quad( x, mybeta ), from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
  }

##------------------------------------
## plot scatterplot of flue vs soilm 
##------------------------------------
# par( las=1, mar=c(4,4,1,1), xaxs="i", yaxs="i" )
  par( las=1 )
  with( 
        mdf,
        heatscatter( 
                    soilm, 
                    flue,
                    ylim=c(0,1.4),
                    xlim=c(0,1),
                    main="",
                    xlab="soil water content (fraction)",
                    ylab="fLUE", 
                    cex.lab=1.2
                    ) 

      )
  box( lwd=1.5 )
  axis( 1, lwd=1.5 ); axis( 1, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
  axis( 2, lwd=1.5, cex.lab=2 ); axis( 2, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
  axis( 3, lwd=1.5, labels=FALSE ); axis( 3, at=seq(0,1.1,by=0.1), tck=-0.02, labels=FALSE )
  axis( 4, lwd=1.5, labels=FALSE ); axis( 4, at=seq(0,1.4,by=0.1), tck=-0.02, labels=FALSE )
  abline( h=1.0, lwd=0.5 )

  for ( sitename in unique(mdf$mysitename) ){
    mybeta <- dplyr::filter( mdf_flue0, mysitename==sitename)$beta
    if (!is.na(mybeta)) mycurve( function(x) stress_quad( x, mybeta ), from=0.0, to=1.0, col='black', add=TRUE, lwd=1 )
  }

  legend( "bottomright", legend=c("low density", "", "", "", "high density"), pch=19, col=colorRampPalette( c("gray65", "navy", "red", "yellow"))(5), bty="n", inset=c(0,0) )
