##------------------------------------------------
## identical to get_linearfit2.R except that two different linear models are fitted - one for grasslands/shrublands and one for others
##------------------------------------------------
get_yintersect <- function( df, target="ratio_obs_mod_pmodel", bin=TRUE, beta_min=0.01, x0_fix=0.8, agg=NA, useweights=FALSE, doplot=FALSE ){

  require(dplyr)
  require(tidyr)
  require(minpack.lm)

  source("stress_quad_1sided.R")

  if (!is.element(target, names(df))||!is.element("soilm_mean", names(df))){
    
    print("ERROR: missing variables fvar or soilm_mean in get_yintersect")
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
    if (useweights){
      weights <- 1.0 - df_tmp$fvar
      weights <- ifelse( weights<0, 0, weights )
      weights <- ifelse( is.na(weights), 0, weights )
      weights <- weights^2
    } else {
      weights <- rep( 1.0, nrow(df_tmp) )
    }

    if (nrow(filter(df_tmp, soilm_mean<0.5))>0){

      ## parabolic stress function
      eq <- paste0( target, " ~ stress_quad_1sided( soilm_mean, x0, beta )")
      fit <- try( 
                  nlsLM( 
                        eq,
                        data=df_tmp,
                        start=list( x0=x0_fix, beta=1.0 ),
                        lower=c( x0_fix,  beta_min ),
                        upper=c( x0_fix,  99  ),
                        algorithm="port",
                        weights = weights
                        )
                  )

      if (doplot){

        par( las=1 )
        with( df_tmp, plot( soilm_mean, fvar, xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("black", 0.2) ) )

        ## Curve
        if (class(fit)!="try-error"){
          ## parabolic stress function
          mycurve(  function(x) stress_quad_1sided( x, x0 = 0.9, coef(fit)[[ "beta" ]] ), from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )
        }
        mtext( df_tmp$mysitename[1], line = 0.5, adj = 0.0, font = 2 )

      }

      # return coefficients of fitted function
      ## parabolic stress function
      if (class(fit)!="try-error"){ 
        out <- c( coef(fit), y0=stress_quad_1sided( 0.0, x0_fix, coef(fit)[[ "beta" ]] ) )
      } else {
        out <- c( x0=NA, beta=NA, y0=NA )
      }
      
    } else {
      out <- c( x0=NA, beta=NA, y0=NA )
    }

    return( out )

  }
}


get_linearfit6 <- function( df, target="ratio_obs_mod_pmodel", monthly=FALSE, bin=TRUE, x0_fix=0.8, agg=NA, useweights=FALSE, doplot=FALSE ){
  ##------------------------------------------------------------------------
  ## This first fits the y-axis intersect using a stress function ('stress_quad_1sided()'),
  ## and then fits a linear model between between this y-axis intersect and the site-level mean alpha.
  ##------------------------------------------------------------------------
  
  require(dplyr)
  require(tidyr)
  require(readr)

  beta_min <- 0.01

  siteinfo <- read_csv("../sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv")

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
                                      dplyr::select( dplyr::filter( df, mysitename==sitename ), mysitename, date, soilm_mean, ratio_obs_mod_pmodel, fvar ), 
                                      target=target, bin=bin, beta_min=beta_min, x0_fix=x0_fix, agg=agg, useweights=useweights, doplot=doplot
                                      ) 
                )
  }

  out <- as.data.frame( cbind( df_flue0, out ) )

  ##------------------------------------------------------------------------
  ## Fit linear models
  ##------------------------------------------------------------------------
  out <- out %>% left_join( select( siteinfo, mysitename, classid ), by = "mysitename" )

  linmod_tree  <- lm( y0 ~ meanalpha, data=dplyr::filter( out, !(classid %in% c("GRA", "CSH") ) & y0 > -1 ) )
  linmod_grass <- lm( y0 ~ meanalpha, data=dplyr::filter( out,   classid %in% c("GRA", "CSH")   & y0 > -1 ) )
  
  return( list( linmod_tree=linmod_tree, linmod_grass=linmod_grass, data=out ) )

}