stress_quad_flattop <- function( x, x0, off, beta ){
  outstress <- 1.0 + off - beta * ( x - x0 ) ^ 2
  outstress <- ifelse( outstress>1, 1, outstress )
  return( outstress )
}

get_yintersect <- function( df, target="ratio_obs_mod", bin=TRUE, beta_min=0.01, x0_fix=0.8 ){

  require( dplyr )
  require(tidyr)
  require(minpack.lm)

  source("calc_flue_est_alpha.R")
  source("calc_flue_est_flue0.R")
  source("stress_quad_1sided.R")

  if (!is.element("fvar", names(df))||!is.element("soilm_mean", names(df))){
    
    print("ERROR: missing variables fvar or soilm_mean in get_yintersect")
    return(NA)

  } else {

    ## Get median by bin and get fvar_vs_soilm for this site (used for clustering)
    if (bin){
      nsbins <- 10
      sbins <- seq( 0.0, 1.0, 1.0/nsbins )
      xvals <- sbins[1:(length(sbins)-1)] + (sbins[2]-sbins[1])/2
      df$insbin <- NULL
      df <- df %>% mutate( insbin = cut( soilm_mean , breaks = sbins ) ) %>% group_by( insbin )
      tmp <- df %>% summarise( median=median( fvar, na.rm=TRUE ) ) %>% complete( insbin, fill = list( median = NA ) ) %>% dplyr::select( median )
      yvals <- unlist(tmp)[1:nsbins]
      df_tmp <- data.frame( soilm_mean=xvals, fvar=yvals )
    } else {
      df_tmp <- df
    }

    # ##------------------------------------------------
    # ## Fit by medians in bis - FLATTOP
    # ##------------------------------------------------
    # fit <- try( 
    #             nlsLM( 
    #                   ratio_obs_mod ~ stress_quad_flattop( soilm_mean, x0, off, beta ),
    #                   data=df_tmp,
    #                   start=list( x0=1.0, off=0.0, beta=1.0 ),
    #                   lower=c( 0.01, -1.0, 0.01 ),
    #                   algorithm="port"
    #                   ) 
    #             )

    ##------------------------------------------------
    ## Fit by medians in bis - 1SIDED
    ##------------------------------------------------
    eq <- paste0( target, "~ stress_quad_1sided( soilm_mean, x0, beta )")
    fit <- try( 
                nlsLM( 
                      # ratio_obs_mod ~ stress_quad_1sided( soilm_mean, x0, beta ),
                      eq,
                      data=df_tmp,
                      start=list( x0=x0_fix, beta=1.0 ),
                      lower=c( x0_fix,  beta_min ),
                      upper=c( x0_fix,  99  ),
                      algorithm="port"
                      ) 
                )

    ## return coefficients of fitted function
    if (class(fit)!="try-error"){ 
      out <- c( coef(fit), y0=stress_quad_1sided( 0.0, x0_fix, coef(fit)[[ "beta" ]] ) )
    } else {
      out <- c( x0=NA, beta=NA, y0=NA )
    }

    return( out )

  }
}


get_linearfit2 <- function( df, target="ratio_obs_mod", monthly=FALSE, bin=TRUE, x0_fix=0.8 ){
  
  require(dplyr)
  require(tidyr)

  beta_min <- 0.01

  if (monthly){
    ## add date and MOY to dataframe nice_agg
    df <- df %>% mutate( date = as.POSIXct( as.Date( paste( as.character( year ), "-01-01", sep="" ) ) + doy - 1 ))
    df <- df %>% mutate( moy = as.numeric( format( date, format="%m" ) ) )

    ## aggregate nice_agg to monthly values
    df <- df %>% group_by( mysitename, year, moy ) %>% summarise( fvar = mean( fvar, na.rm=TRUE ), soilm_mean = mean( soilm_mean, na.rm=TRUE ) )    
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
    # df_flue0$flue0[ which(df_flue0$mysitename==sitename) ] <- get_yintersect( dplyr::filter( df, mysitename==sitename ) )
    out <- rbind( out, get_yintersect( dplyr::filter( df, mysitename==sitename ), target=target, bin=bin, beta_min=beta_min, x0_fix=x0_fix ) )
  }

  out <- as.data.frame( cbind( df_flue0, out ) )

  ##------------------------------------------------------------------------
  ## Fit linear model
  ##------------------------------------------------------------------------
  linmod <- lm( y0 ~ meanalpha, data=dplyr::filter( out, y0 > -1 ) )
  # linmod <- lm( y0 ~ meanalpha, data=out )
  
  return( list( linmod=linmod, data=out ) )

}