get_yintersect <- function( df, target="ratio_obs_mod_pmodel", bin=TRUE, beta_min=0.01, x0_fix=0.8, agg=NA ){

  require(dplyr)
  require(tidyr)
  require(minpack.lm)

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
    eq <- paste0( target, " ~ stress_quad_1sided( soilm_mean, x0, beta )")
    fit <- try( 
                nlsLM( 
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


get_linearfit2 <- function( df, target="ratio_obs_mod_pmodel", monthly=FALSE, bin=TRUE, x0_fix=0.8, agg=NA ){
  
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
    out <- rbind( out, get_yintersect( 
                                      dplyr::select( dplyr::filter( df, mysitename==sitename ), year, doy, soilm_mean, ratio_obs_mod_pmodel, fvar ), 
                                      target=target, bin=bin, beta_min=beta_min, x0_fix=x0_fix, agg=agg 
                                      ) 
                )
  }

  out <- as.data.frame( cbind( df_flue0, out ) )

  ##------------------------------------------------------------------------
  ## Fit linear model
  ##------------------------------------------------------------------------
  linmod <- lm( y0 ~ meanalpha, data=dplyr::filter( out, y0 > -1 ) )
  
  return( list( linmod=linmod, data=out ) )

}