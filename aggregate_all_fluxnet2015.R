.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

library( dplyr )
library( cgwtools )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/utils_sofun/analysis_sofun/remove_outliers.R", sep="" ) )
source( "./func_flue_est.R" )

## IMPORTANT: USE SOILMOISTURE FROM S13 FOR NN-TRAINING
# load( paste( myhome, "data/fluxnet_sofun/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata", sep="" ) ) # "new data" with s13

# func_test <- function( x, add ){ 
#   y <- x+0.1+add
#   return(y) 
# }


##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
siteinfo <- read.csv( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/soilm_data_usability_fluxnet2015.csv", sep="" ), as.is=TRUE )
do.sites <- dplyr::filter( siteinfo, code!=0 )$mysitename
nam_target  = "lue_obs_evi"
use_weights = FALSE
use_fapar   = FALSE


## Manual settings ----------------
# do.sites   = "AU-Dry"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
overwrite_nice = TRUE
overwrite_modis = TRUE
overwrite_mte = TRUE
verbose    = FALSE
##---------------------------------


##------------------------------------------------
## Get MTE-GPP for all sites
##------------------------------------------------
filn <- paste( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_8Days_4Beni.csv", sep="" )
if ( file.exists( filn ) ){
  mte_dl <- read.csv( paste( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_Daily_4Beni.csv", sep="" ), as.is=TRUE )
  mte_8d <- read.csv( filn, as.is=TRUE )

  ## replace with NA
  for (ivar in names(mte_8d)){ mte_8d[[ ivar ]][ which( mte_8d[[ ivar ]]==-9999 ) ] <- NA }
  # for (ivar in names(dmte)){
  #   dmte[[ ivar ]][ which( dmte[[ ivar ]]==-9999 ) ] <- NA
  # }

  avl_data_mte <- TRUE
} else {
  avl_data_mte <- FALSE
}

##------------------------------------------------
## Check availability of MODIS GPP data
##------------------------------------------------
fillist <- list.files( paste( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/", sep="" ) )
if (length(fillist)==0) {
  avl_data_modis <- FALSE
} else {
  avl_data_modis <- TRUE
}


## check and override if necessary
if ( nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
  plotlue <- TRUE
  if (nam_target=="lue_obs_evi"){
    fapar_data <- "evi"
  } else if (nam_target=="lue_obs_fpar"){
    fapar_data <- "fpar"
  }
  if (use_fapar){
    print("WARNING: setting use_fapar to FALSE")
    use_fapar <- FALSE
  }
}

## identifier for output files
if (use_fapar){
  if (nam_target=="lue_obs_evi"){
    char_fapar <- "_withEVI"
  } else if (nam_target=="lue_obs_fpar"){
    char_fapar <- "_withFPAR"
  } else {
    print("ERROR: PROVIDE VALID FAPAR DATA!")
  }
} else {
  char_fapar <- ""
}

if (use_weights){
  char_wgt <- "_wgt"
} else {
  char_wgt <- ""
}

print( "Aggregating and complementing data for all sites ..." )

##------------------------------------------------
## Initialise aggregated data
##------------------------------------------------
## fvar and soilm data to be complemented with cluster info
nice_agg  <- data.frame()
mte_agg   <- data.frame()
modis_agg <- data.frame()
# nice_resh         <- data.frame()

## all possible soil moisture datasets
varnams_swc_full <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob", "soilm_obs" )

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  nicefiln <- paste0("data/nice_all_", sitename, ".Rdata" )

  if (file.exists(nicefiln)&&!overwrite_nice){

    load( nicefiln )

  } else {

    ##------------------------------------------------
    ## load site data and "detatch"
    ##------------------------------------------------
    nice <- fluxnet[[ sitename ]]$ddf$s13

    if (!is.null(nice)){

      nice <- nice %>% dplyr::select( year_dec, year, doy, moy, dom, soilm_splash150=wcont, gpp_pmodel=gpp, aet_pmodel=aet, pet_pmodel=pet )

      ##------------------------------------------------
      ## Get alternative soil moisture data
      ##------------------------------------------------
      nice$soilm_splash220 <- fluxnet[[ sitename ]]$ddf$s11$wcont
      nice$soilm_swbm      <- fluxnet[[ sitename ]]$ddf$s12$wcont
      nice$soilm_etobs     <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et
      nice$soilm_etobs_ob  <- fluxnet[[ sitename ]]$ddf$swc_by_etobs$soilm_from_et_orthbucket

      varnams_swc <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob" )

      ##------------------------------------------------
      ## normalise soil moisture
      ##------------------------------------------------
      # nice <- nice %>% mutate( soilm_splash150 = soilm_splash150 / 150 )
      # nice <- nice %>% mutate( soilm_splash220 = soilm_splash220 / 220 )
      # nice <- nice %>% mutate( soilm_swbm      = soilm_swbm      / 220 )
      # nice <- nice %>% mutate( soilm_etobs     = soilm_etobs     / 220 )
      # nice <- nice %>% mutate( soilm_etobs_ob  = soilm_etobs_ob  / 220 )

      nice <- nice %>% mutate( soilm_splash150 = soilm_splash150 / max( soilm_splash150, na.rm=TRUE ) )
      nice <- nice %>% mutate( soilm_splash220 = soilm_splash220 / max( soilm_splash220, na.rm=TRUE ) )
      nice <- nice %>% mutate( soilm_swbm      = soilm_swbm      / max( soilm_swbm     , na.rm=TRUE ) )
      nice <- nice %>% mutate( soilm_etobs     = soilm_etobs     / max( soilm_etobs    , na.rm=TRUE ) )
      nice <- nice %>% mutate( soilm_etobs_ob  = soilm_etobs_ob  / max( soilm_etobs_ob , na.rm=TRUE ) )

      ##------------------------------------------------
      ## Get observational soil moisture data (variable availability!)
      ##------------------------------------------------
      varnams_swc_obs <- c()

      ## Get observational soil moisture data (variable availability!)
      if ( is.element( "soilm_obs", varnams_swc ) ){
        relevant <- names(fluxnet[[ sitename ]]$ddf$swc_obs)[(is.element( substr(names(fluxnet[[ sitename ]]$ddf$swc_obs), start=1, stop=3), "SWC" ))]
        if (length(relevant)>0){
          for (iobs in seq(length(relevant))){
            varnam <- paste( "soilm_obs_", iobs, sep="" )
            nice[[ varnam ]] <- fluxnet[[ sitename ]]$ddf$swc_obs[[ relevant[iobs] ]]
            varnams_swc_obs <- c( varnams_swc_obs, varnam )
          }
          varnams_swc <- c( varnams_swc, "soilm_obs" )
        } else {
          ## no obs. soilm. data found, remove from data list
          varnams_swc <- varnams_swc[ -which( "soilm_obs" ) ]
        }  
      }

    }

    ## Average over soil moisture datasets; all and obs-only
    if (length(varnams_swc_obs)>0){
      nice$soilm_obs <- apply( dplyr::select( nice, one_of(varnams_swc_obs)), 1, FUN=mean, na.rm=TRUE )
      nice$soilm_obs[ is.nan( nice$soilm_obs ) ] <- NA
    }

    nice$soilm_mean <- apply( dplyr::select( nice, one_of(varnams_swc)), 1, FUN=mean, na.rm=TRUE )
    nice$soilm_mean[ is.nan( nice$soilm_mean ) ] <- NA

    ##------------------------------------------------
    ## Get observational data and add to 'data'
    ##------------------------------------------------
    obs <- dplyr::select( fluxnet[[ sitename ]]$ddf$obs, year, moy, dom, gpp_obs2015_GPP_NT_VUT_REF, gpp_obs2015_GPP_NT_VUT_REF_gfd, le_f_mds )

    nice <- nice  %>% left_join( obs, by=c( "year", "moy", "dom" ) ) %>%
                      rename( gpp_obs=gpp_obs2015_GPP_NT_VUT_REF, gpp_obs_gfd=gpp_obs2015_GPP_NT_VUT_REF_gfd, et_obs=le_f_mds ) %>%
                      mutate( wue_obs=gpp_obs/(et_obs*1e-6) )

    ##------------------------------------------------
    ## Get input data
    ##------------------------------------------------
    inp <- fluxnet[[ sitename ]]$ddf$inp %>% dplyr::select( year, moy, dom, temp, ppfd, vpd, prec, evi, fpar )

    nice <- nice %>% left_join( inp, by=c( "year", "moy", "dom" ) )

    ##------------------------------------------------
    ## get LUE and remove outliers
    ##------------------------------------------------
    nice <- nice %>% mutate( lue_obs_evi  = remove_outliers( gpp_obs / ( ppfd * evi  ), coef=3.0 ) )
    nice <- nice %>% mutate( lue_obs_fpar = remove_outliers( gpp_obs / ( ppfd * fpar ), coef=3.0 ) )      

    for (ivar in varnams_swc_obs){
      nice[[ ivar ]]  <- nice[[ ivar ]]  / max( nice[[ ivar ]] , na.rm=TRUE )
    }

    ##------------------------------------------------
    ## additional variables
    ##------------------------------------------------
    nice$bias_pmodel  <-  nice$gpp_pmodel / nice$gpp_obs
    nice$bias_pmodel[ which(is.infinite(nice$bias_pmodel)) ] <- NA

    nice$ratio_obs_mod  <-  nice$gpp_obs / nice$gpp_pmodel
    nice$ratio_obs_mod[ which(is.infinite(nice$ratio_obs_mod)) ] <- NA

    nice$alpha  <-  nice$aet_pmodel / nice$pet_pmodel
    nice$alpha[ which(is.infinite(nice$alpha)) ] <- NA

    ## add row to aggregated data
    nice <- nice %>% mutate( mysitename=sitename )

    ## fLUE estimate based on current soil moisture and average AET/PET
    meanalpha <- mean( nice$aet_pmodel / nice$pet_pmodel, na.rm=TRUE )
    nice <- nice %>%  mutate( flue0    = calc_flue0( meanalpha ) )%>% 
                      mutate( beta = calc_beta( flue0 ) ) %>% 
                      mutate( flue_est = func_flue_est( soilm_mean, beta ) ) %>% 
                      mutate( dry=ifelse(alpha<1, TRUE, FALSE))

    ##------------------------------------------------
    ## save to file
    ##------------------------------------------------
    save( nice, file=nicefiln )

  }

  ##------------------------------------------------
  ## record for aggregated data
  ##------------------------------------------------
  usecols <- c(
                "year_dec",
                "year",
                "doy",
                "gpp_obs",
                "gpp_pmodel",
                "aet_pmodel",
                "pet_pmodel",
                "alpha",
                "ppfd",
                "vpd",
                "evi", 
                "fpar", 
                "soilm_splash220",
                "soilm_splash150",
                "soilm_swbm",
                "soilm_mean",
                "bias_pmodel", 
                "ratio_obs_mod", 
                "lue_obs_evi", 
                "lue_obs_fpar",
                "flue_est",
                "dry"
                )

  # nice_agg <- dplyr::select( nice, one_of( usecols ) ) %>% cbind( mysitename, . ) %>% rbind( . )

  nice_agg <- rbind( nice_agg, dplyr::select( nice, one_of( usecols ) ) )

  # ##------------------------------------------------
  # ## Reshape dataframe to stack data from differen soil moisture datasets along rows
  # ##------------------------------------------------
  # for (isoilm in varnams_swc_full){
  #   if (isoilm=="soilm_obs") nice$soilm_obs <- nice$soilm_obs_mean
  #   if (is.null(nice[[ isoilm ]])){
  #     nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ paste( "moist_",      isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
  #     nice[[ isoilm ]]                                 <- rep( NA, nrow(nice) )
  #   }
  #   addrows <- data.frame( 
  #                           mysitename = sitename,
  #                           year_dec   = nice$year_dec,
  #                           var_nn_act = nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]],
  #                           var_nn_pot = nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]],
  #                           var_nn_vpd = nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]],
  #                           moist      = nice[[ paste( "moist_",      isoilm, sep="" ) ]],
  #                           lue_obs_evi= nice$lue_obs_evi,
  #                           soilm      = nice[[ isoilm ]],
  #                           vpd        = nice$vpd,
  #                           iabs       = nice$evi * nice$ppfd,
  #                           soilm_data = rep( isoilm, nrow(nice) )
  #                           )
  #   nice_resh <- rbind( nice_resh, addrows )
  # }


  ##------------------------------------------------
  ## Get MTE-GPP for this site
  ##------------------------------------------------
  if (avl_data_mte){

    if (is.element( sitename, mte_8d$Site.code)){

      filn <- paste( "data/mte_", sitename, ".Rdata", sep="" )

      if ( file.exists(filn) && !overwrite_mte ){

        ## load 'nice_to_mte' from file
        load( filn )

      } else {

        ## prepare dataframe 'nice_to_mte'
        # dmte <- dplyr::filter( mte_dl, Site.code==sitename )
        missing_mte <- FALSE

        ## make mte a bit nicer
        mte <-  filter( mte_8d, Site.code==sitename ) %>%
                rename( mysitename=Site.code, gpp_mte=MTE, gpp_mte_m=MTE_M, gpp_mte_viterbo=MTE_Viterbo, gpp_rf=RF, gpp_rf_fromdaily=RF_from_daily, doy_start=StartDoY, doy_end=EndDoY, year_start=StartYear, year_end=EndYear ) %>%
                mutate( doy_end = doy_end - 1 ) %>%  # assuming that the doy_end is not counted towards bins aggregate
                mutate( date_start = as.POSIXct( as.Date( paste( as.character(year_start), "-01-01", sep="" ) ) + doy_start - 1 ),
                        date_end   = as.POSIXct( as.Date( paste( as.character(year_end  ), "-01-01", sep="" ) ) + doy_end   - 1 ) )

        ## group by 8d bins from MTE data
        nice <- nice %>%  mutate( date = as.POSIXct( as.Date( paste( as.character(year), "-01-01", sep="" ) ) + doy - 1 ) ) %>% 
                          # mutate( inmtebin = cut( as.numeric(date), breaks = c( mte$date_start[1], mte$date_end ), include.lowest = TRUE ) ) 
                          mutate( inmtebin = cut( as.numeric(date), breaks = c( mte$date_start ), right = FALSE ) ) 

        ## summarise by bin taking means
        nice_to_mte <- nice %>% group_by( inmtebin ) %>% summarise_all( mean, na.rm=TRUE )
        tmp         <- nice %>% select( inmtebin, doy, year ) %>% group_by( inmtebin ) %>% summarise( doy_start=min( doy ), doy_end=max( doy ), year_start=min( year ), year_end=max( year ) )
        nice_to_mte <- cbind( nice_to_mte, select( tmp, doy_start, doy_end, year_start, year_end ) )

        ## merge dataframes (averaged nice and mte)
        nice_to_mte <- nice_to_mte %>% select( doy_start, year_start, one_of( usecols ) ) %>% left_join( mte, by=c("doy_start", "year_start") )
        
        ## get additional variables
        nice_to_mte <- nice_to_mte %>%  mutate( bias_mte = gpp_mte / gpp_obs )          %>% mutate( bias_mte=ifelse( is.infinite( bias_mte ), NA, bias_mte ) ) %>% 
                                        mutate( ratio_obs_mod_mte = gpp_obs / gpp_mte ) %>% mutate( ratio_obs_mod_mte=ifelse( is.infinite( bias_mte ), NA, ratio_obs_mod_mte ) ) %>%  
                                        mutate( bias_rf = gpp_rf / gpp_obs )            %>% mutate( bias_rf=ifelse( is.infinite( bias_rf ), NA, bias_rf ) )  %>%  
                                        mutate( ratio_obs_mod_rf = gpp_obs / gpp_rf )   %>% mutate( ratio_obs_mod_rf=ifelse( is.infinite( bias_rf ), NA, ratio_obs_mod_rf ) )

        ## save to file
        save( nice_to_mte, file=filn )

      }

      ## add row to aggregated data
      mte_agg <- rbind( mte_agg, nice_to_mte )

    } else {

      missing_mte <- TRUE

    }

  }


  ##------------------------------------------------
  ## Get MODIS-GPP for this site
  ##------------------------------------------------
  if (avl_data_modis){

    filn <- paste( myhome, "data/modis_", sitename, ".Rdata", sep="" )

    if ( file.exists(filn) && !overwrite_modis ){

      ## load 'nice_to_modis' from file
      load( filn )
      avl_modisgpp <- TRUE

    } else {

      ## prepare 'nice_to_modis'
      modis <- try( read.csv( paste( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/", sitename, "/gpp_8d_modissubset_", sitename, ".csv", sep="" ), as.is=TRUE ))
      if (class(modis)!="try-error"){
        avl_modisgpp <- TRUE

        ## make modis a bit nicer
        modis <- modis %>% rename( gpp_modis = data ) %>% mutate( gpp_modis = gpp_modis / 8.0, doy_start = doy, doy_end = lead( doy ) - 1, year_start = year, date_start = as.POSIXct( as.Date( date ) ), date_end = as.POSIXct( as.Date( lead( date ) ) - 1 ) )

        ## group nice by 8d bins from MODIS data
        nice <- nice %>%  mutate( date = as.POSIXct( as.Date( paste( as.character(year), "-01-01", sep="" ) ) + doy - 1 ) ) %>% 
                          mutate( inmodisbin = cut( as.numeric(date), breaks = c( modis$date_start ), right = FALSE ) )

        ## summarise by bin taking means
        nice_to_modis <- nice %>% group_by( inmodisbin ) %>% summarise_all( mean, na.rm=TRUE )
        tmp           <- nice %>% select( inmodisbin, doy, year ) %>% group_by( inmodisbin ) %>% summarise( doy_start=min( doy ), doy_end=max( doy ), year_start=min( year ), year_end=max( year ) )
        nice_to_modis <- cbind( nice_to_modis, select( tmp, doy_start, doy_end, year_start, year_end ) )

        ## merge dataframes (averaged nice and modis)
        nice_to_modis <- nice_to_modis %>% select( doy_start, year_start, one_of( usecols ) ) %>% left_join( modis, by=c("doy_start", "year_start") )
        
        ## get additional variables
        nice_to_modis <- nice_to_modis %>% mutate( bias_modis = gpp_modis / gpp_obs )          %>% mutate( bias_modis=ifelse( is.infinite( bias_modis ), NA, bias_modis ) ) %>% 
                                           mutate( ratio_obs_mod_modis = gpp_obs / gpp_modis ) %>% mutate( ratio_obs_mod_modis=ifelse( is.infinite( bias_modis ), NA, ratio_obs_mod_modis ) )

        ## save to file
        save( nice_to_modis, file=filn )

      } else {

        avl_modisgpp <- FALSE

      }

    }

    ## add row to aggregated data
    if (avl_modisgpp){
      modis_agg <- rbind( modis_agg, nice_to_modis )
    }

  }

}

print("... done.")

if ( length( dplyr::filter( siteinfo, code!=0 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  save( nice_agg,  file=paste("data/nice_all_agg_",  nam_target, char_fapar, ".Rdata", sep="") )
  # save( nice_all_resh, file=paste("data/nice_all_resh_", nam_target, char_fapar, ".Rdata", sep="") )

  if (avl_data_mte)   save( nice_to_mte_agg,   file=paste("data/nice_all_mte_agg_",   nam_target, char_fapar, ".Rdata", sep="") )
  if (avl_data_modis) save( nice_to_modis_agg, file=paste("data/nice_all_modis_agg_", nam_target, char_fapar, ".Rdata", sep="") )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

