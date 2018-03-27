## This works like get_linearfit5.R but with a different functional form of the soil moisture stress function: 2-parameter exponential instead of parabolic

get_yintersect <- function( df, target="fvar", bin=TRUE, beta_min=0.01, x0_fix=0.8, agg=NA, useweights=FALSE, doplot=FALSE ){

  require(dplyr)
  require(tidyr)
  require(minpack.lm)

  source("stress_exp.R")

  if (!is.element(target, names(df))||!is.element("soilm_splash", names(df))){
    
    print("ERROR: missing variables fvar or soilm_splash in get_yintersect")
    return(NA)

  } else {

    ## Get median by bin and get fvar_vs_soilm for this site (used for clustering)
    if (bin){
      nbins <- 10
      bins <- seq( 0.0, 1.0, 1.0/nbins )
      xvals <- bins[1:(length(bins)-1)] + (bins[2]-bins[1])/2
      df$inbin <- NULL
      df <- df %>% mutate( inbin = cut( fvar , breaks = bins ) ) %>% group_by( inbin )
      tmp <- df %>% summarise( median=median( fvar, na.rm=TRUE ) ) %>% complete( inbin, fill = list( median = NA ) ) %>% dplyr::select( median )
      yvals <- unlist(tmp)[1:nbins]
      df_tmp <- data.frame( soilm_splash=xvals, fvar=yvals )
    } else {
      df_tmp <- df
    }

    ##------------------------------------------------
    ## Aggregate data to N-daily means before fitting
    ##------------------------------------------------
    if (!is.na(agg)){
      df_tmp <- df_tmp %>%  mutate( date = as.POSIXct( as.Date( paste( as.character(year), "-01-01", sep="" ) ) + doy - 1 ) )
      breaks <- seq.POSIXt( df_tmp$date[1], df_tmp$date[nrow(df_tmp)], by=paste0( as.character(agg), " days" ) )
      df_tmp <- df_tmp %>%  mutate( inbin = cut( as.numeric(date), breaks = breaks ) ) %>%  #, right = FALSE
                            group_by( inbin ) %>% summarise_all( mean, na.rm=TRUE )
    }

    ##------------------------------------------------
    ## Fit by medians in bis - 1SIDED
    ##------------------------------------------------
    if (useweights){
      weights <- 1.0 - df_tmp$fvar
      weights <- ifelse( weights<0, 0, weights )
      weights <- ifelse( is.na(weights), 0, weights )
      weights <- weights^2
    } else {
      weights <- rep( 1.0, nrow(df_tmp) )
    }

    ## exponential stress function
    eq <- paste0( target, " ~ stress_exp( soilm_splash, y0, curve )")
    fit <- try( 
                nlsLM( 
                      eq,
                      data=df_tmp,
                      start=list( y0=0.0, curve=3.0 ),
                      lower=c( -5,  1.0 ),
                      upper=c( 1.0, 99.0  ),
                      algorithm="port"
                      )
                )

    if (doplot){

      par( las=1 )
      with( df_tmp, plot( soilm_splash, fvar, xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("black", 0.2) ) )

      ## Curve
      if (class(fit)!="try-error"){
        ## exponential stress function
        mycurve(  function(x) stress_exp( x, y0 = coef(fit)[[ "y0" ]], curve = coef(fit)[[ "curve" ]]), from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )
      }
      mtext( df_tmp$mysitename[1], line = 0.5, adj = 0.0, font = 2 )

    }
    
    ## exponential stress function
    if (class(fit)!="try-error"){ 
      out <- c( coef(fit) )
    } else {
      out <- c( y0=NA, curve=NA )
    }    

    return( out )

  }
}


get_linearfit5 <- function( df, target="ratio_obs_mod_pmodel", monthly=FALSE, bin=TRUE, x0_fix=0.8, agg=NA, useweights=FALSE, doplot=FALSE ){
  ##------------------------------------------------------------------------
  ## This first fits the y-axis intersect using a stress function ('stress_quad_1sided()'),
  ## and then fits a linear model between between this y-axis intersect and the site-level mean alpha.
  ##------------------------------------------------------------------------
  
  require(dplyr)
  require(tidyr)

  beta_min <- 0.01

  if (monthly){
    ## add date and MOY to dataframe nice_agg
    df <- df %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
    df <- df %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

    ## aggregate nice_agg to monthly values
    df <- df %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_splash = mean( soilm_splash, na.rm=TRUE ) )    
  }

  df_flue0 <- df %>% dplyr::select( mysitename ) %>% unique()

  ## Merge mean annual alpha (AET/PET) values into this dataframe
  load( "../sofun/utils_sofun/analysis_sofun/fluxnet2015/data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
  df_flue0 <- df_flue0 %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

  ##------------------------------------------------------------------------
  ## Get y-axis intersect for each site
  ##------------------------------------------------------------------------
  df_flue0 <- df_flue0 %>% mutate( flue0 = NA )
  out <- c()
  for (sitename in df_flue0$mysitename){
    out <- rbind( out, get_yintersect( 
                                      dplyr::select( dplyr::filter( df, mysitename==sitename ), mysitename, date, soilm_splash, ratio_obs_mod_pmodel, fvar ), 
                                      target=target, bin=bin, beta_min=beta_min, x0_fix=x0_fix, agg=agg, useweights=useweights, doplot=doplot
                                      ) 
                )
  }

  out <- as.data.frame( cbind( df_flue0, out ) )

  ##------------------------------------------------------------------------
  ## Fit linear models
  ##------------------------------------------------------------------------
  linmod_y0 <- lm( y0 ~ meanalpha, data=out )
  # linmod_curve <- lm( curve ~ y0, data=out )
  linmod_curve <- lm( curve ~ meanalpha, data=out )
  
  return( list( linmod_curve=linmod_curve, linmod_y0=linmod_y0, data=out ) )

}