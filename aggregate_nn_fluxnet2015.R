.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

library( dplyr )
library( cgwtools )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/utils_sofun/analysis_sofun/remove_outliers.R", sep="" ) )

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/successcodes.csv", sep="" ), as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
# do.sites   = "AU-Rob"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
overwrite_modis = FALSE
overwrite_mte = FALSE
verbose    = FALSE
##---------------------------------

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

##------------------------------------------------
## Get MTE-GPP for all sites
##------------------------------------------------
filn <- paste( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_8Days_4Beni.csv", sep="" )
if ( file.exists( filn ) ){
  mte_dl <- read.csv( paste( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_Daily_4Beni.csv", sep="" ), as.is=TRUE )
  mte_8d <- read.csv( filn, as.is=TRUE )
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
nice_agg          <- data.frame()
nice_to_mte_agg   <- data.frame()
nice_to_modis_agg <- data.frame()
nice_resh         <- data.frame()

## all possible soil moisture datasets
varnams_swc_full <- c( "soilm_splash150", "soilm_splash220", "soilm_swbm", "soilm_etobs", "soilm_etobs_ob", "soilm_obs" )

jdx <- 0
for (sitename in do.sites){

  jdx <- jdx + 1
  missing_mte <- FALSE

  infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 
  if (verbose) print( paste( "opening file", infil ) )

  ##------------------------------------------------
  ## load nn_fVAR data and "detatch"
  ##------------------------------------------------
    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

    nice$bias_pmodel  <-  nice$gpp_pmodel / nice$gpp_obs
    nice$bias_pmodel[ which(is.infinite(nice$bias_pmodel)) ] <- NA

    nice$ratio_obs_mod  <-  nice$gpp_obs / nice$gpp_pmodel
    nice$ratio_obs_mod[ which(is.infinite(nice$ratio_obs_mod)) ] <- NA

    ## add row to aggregated data
    mysitename <- data.frame( mysitename=rep( sitename, nrow(nice) ) )

    if (nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar"){
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act * evi * ppfd, gpp_nn_pot = var_nn_pot * evi * ppfd, gpp_nn_vpd = var_nn_vpd * evi * ppfd )
    } else {
      nice <- nice %>% mutate( gpp_nn_act = var_nn_act, gpp_nn_pot = var_nn_pot, gpp_nn_vpd = var_nn_vpd )
    }

    ## fill with NA if respective soil moisture data was not used
    for (isoilm in varnams_swc_full){
      if ( !is.element( paste( "var_nn_act_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_act_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "var_nn_pot_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_pot_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "var_nn_vpd_", isoilm, sep="" ), names(nice) ) ) nice[[ paste( "var_nn_vpd_", isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
      if ( !is.element( paste( "moist_",      isoilm, sep="" ), names(nice) ) ) nice[[ paste( "moist_",      isoilm, sep="" ) ]] <- rep( NA, nrow(nice) )
    }


  ##------------------------------------------------
  ## Re-save data after additions to 'nice' dataframe
  ##------------------------------------------------
  # if (verbose) print( paste( "resaving nice into file", infil ) )
  # if (verbose) print( "names:" ); print( names(nice) )
  nn_fluxnet[[ sitename ]]$nice <- nice
  resave( nn_fluxnet, file=infil )

  ##------------------------------------------------
  ## record for aggregated data
  ##------------------------------------------------
  usecols <- c(
                "year_dec",
                "year",
                "doy",
                "gpp_obs",
                "var_nn_pot", 
                "var_nn_act",
                "var_nn_vpd",
                "gpp_nn_act", 
                "gpp_nn_pot",
                "fvar", 
                "gpp_pmodel",
                "aet_pmodel",
                "pet_pmodel",
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
                "is_drought_byvar", 
                "lue_obs_evi", 
                "lue_obs_fpar",
                paste( "var_nn_act_", varnams_swc_full, sep="" ),
                paste( "var_nn_pot_", varnams_swc_full, sep="" ),
                paste( "var_nn_vpd_", varnams_swc_full, sep="" ),
                paste( "moist_", varnams_swc_full, sep="" )
                )

  # nice_agg <- dplyr::select( nice, one_of( usecols ) ) %>% cbind( mysitename, . ) %>% rbind( . )

  sub <- dplyr::select( nice, one_of( usecols ) )
  nice_agg <- rbind( nice_agg, cbind(  mysitename, sub ) )


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


  # ##------------------------------------------------
  # ## Get MTE-GPP for this site
  # ##------------------------------------------------
  # if (avl_data_mte){

  #   if (is.element( sitename, mte_8d$Site.code)){

  #     filn <- paste( myhome, "data/mte_", sitename, ".Rdata", sep="" )

  #     if ( file.exists(filn) && !overwrite_mte ){

  #       ## load 'nice_to_mte' from file
  #       load( filn )

  #     } else {

  #       ## prepare dataframe 'nice_to_mte'
  #       mte  <- dplyr::filter( mte_8d, Site.code==sitename )
  #       dmte <- dplyr::filter( mte_dl, Site.code==sitename )
  #       missing_mte <- FALSE

  #       mte$year_dec_start <- mte$StartYear + ( mte$StartDoY - 1 ) / 365 
  #       mte$year_dec_end   <- mte$EndYear   + ( mte$EndDoY   - 1 ) / 365 
  #       mte$year_dec       <- mte$StartYear + ( (mte$StartDoY + mte$EndDoY)/2 - 1 ) / 365 ## year_dec is the start of the period. approx( ..., method="linear") holds value [i] constant from year_dec[1] to year_dec[i+1]

  #       dmte$year_dec <- dmte$StartYear + ( dmte$StartDoY - 1 ) / 365 

  #       mte  <- dplyr::rename( mte,  mysitename=Site.code, gpp_mte=MTE, gpp_mte_m=MTE_M, gpp_mte_viterbo=MTE_Viterbo, gpp_rf=RF, gpp_rf_fromdaily=RF_from_daily )
  #       dmte <- dplyr::rename( dmte, mysitename=Site.code, doy=StartDoY, year=StartYear, gpp_rf=RFdaily )
  #       # mte <- dplyr::rename( mte, mysitename=Site.code, gpp_mte=MTE )

  #       ## replace with NA
  #       for (ivar in names(mte)){
  #         mte[[ ivar ]][ which( mte[[ ivar ]]==-9999 ) ] <- NA
  #       }
  #       for (ivar in names(dmte)){
  #         dmte[[ ivar ]][ which( dmte[[ ivar ]]==-9999 ) ] <- NA
  #       }


  #       ## Make 'nice' dataframe conform with 'mte'
  #       nice_to_mte <- c()
  #       for (idx in 1:nrow(mte)){
  #         sub <- dplyr::filter( nice, year_dec>=mte$year_dec_start[idx] & year_dec<=mte$year_dec_end[idx] )
  #         year_dec_save <- sub$year_dec[1]
  #         addline <- unlist( unname( apply( sub, 2, FUN=mean, na.rm=TRUE )))
  #         # addline[1] <- year_dec_save
  #         nice_to_mte <- rbind( nice_to_mte, addline )
  #       }

  #       nice_to_mte <- as.data.frame( nice_to_mte )
  #       colnames( nice_to_mte ) <- names( sub )
  #       rownames( nice_to_mte ) <- NULL

  #       mycolnames <- c( "gpp_mte", "gpp_mte_m", "gpp_mte_viterbo", "gpp_rf", "gpp_rf_fromdaily" )

  #       for (ivar in mycolnames){
  #         nice_to_mte[[ ivar ]] <- mte[[ ivar ]]
  #       }

  #       nice_to_mte$is_drought_byvar <- with( nice_to_mte, ifelse( is_drought_byvar<0.5, FALSE, TRUE ) )

  #       ## get bias measure: ( mod / obs )
  #       ## MTE
  #       nice_to_mte$bias_mte <- nice_to_mte$gpp_mte / nice_to_mte$gpp_obs
  #       nice_to_mte$bias_mte[ which( is.infinite( nice_to_mte$bias_mte ) ) ] <- NA

  #       nice_to_mte$ratio_obs_mod_mte  <-  nice_to_mte$gpp_obs / nice_to_mte$gpp_mte
  #       nice_to_mte$ratio_obs_mod_mte[ which(is.infinite(nice_to_mte$ratio_obs_mod_mte)) ] <- NA        

  #       ## RF
  #       nice_to_mte$bias_rf <- nice_to_mte$gpp_rf / nice_to_mte$gpp_obs
  #       nice_to_mte$bias_rf[ which( is.infinite( nice_to_mte$bias_rf ) ) ] <- NA

  #       nice_to_mte$ratio_obs_mod_rf  <-  nice_to_mte$gpp_obs / nice_to_mte$gpp_rf
  #       nice_to_mte$ratio_obs_mod_rf[ which(is.infinite(nice_to_mte$ratio_obs_mod_rf)) ] <- NA        


  #       ## add dmte to 'nice' dataframe
  #       nice$gpp_rf_daily <- rep( NA, nrow(nice) )
  #       for (idx in 1:nrow(nice)){
  #         idx_use <- which.min( abs( nice$year_dec[idx] - dmte$year_dec ) )
  #         nice$gpp_rf_daily[idx] <- dmte$gpp_rf[ idx_use ]

  #         ## add tolerance
  #         if ( (nice$year_dec[idx] - dmte$year_dec[idx_use]) > 1.5/365 ) nice$gpp_rf_daily[idx] <- NA
  #       }

  #       ## get bias measure: ( mod / obs )
  #       nice$bias_rf <- nice$gpp_rf_daily / nice$gpp_obs
  #       nice$bias_rf[ which( is.infinite( nice$bias_rf ) ) ] <- NA

  #       nice$ratio_obs_mod_dmte  <-  nice$gpp_obs / nice$gpp_rf_daily
  #       nice$ratio_obs_mod_dmte[ which(is.infinite(nice$ratio_obs_mod_dmte)) ] <- NA        


  #       # ## Make 'nice' dataframe conform with 'mte'
  #       # nice_to_dmte <- c()
  #       # for (idx in 1:nrow(mte)){
  #       #   sub <- dplyr::filter( nice, year_dec==mte$year_dec[idx] )
  #       #   year_dec_save <- sub$year_dec[1]
  #       #   addline <- unlist( unname( apply( sub, 2, FUN=mean, na.rm=TRUE ) ) )
  #       #   # addline[1] <- year_dec_save
  #       #   nice_to_dmte <- rbind( nice_to_dmte, addline )
  #       # }

  #       # nice_to_dmte <- as.data.frame( nice_to_dmte )
  #       # colnames( nice_to_dmte ) <- names( sub )
  #       # rownames( nice_to_dmte ) <- NULL

  #       # mycolnames <- c( "gpp_rf" )

  #       # for (ivar in mycolnames){
  #       #   nice_to_dmte[[ ivar ]] <- mte[[ ivar ]]
  #       # }

  #       # nice_to_dmte$is_drought_byvar <- with( nice_to_dmte, ifelse( is_drought_byvar<0.5, FALSE, TRUE ) )

  #       # ## get bias measure: log( mod / obs )
  #       # nice_to_dmte$bias_mte <- nice_to_dmte[[ ivar ]] / nice_to_dmte$gpp_obs
  #       # nice_to_dmte$bias_mte[ which( is.infinite( nice_to_dmte$bias_mte ) ) ] <- NA

  #       # nice_to_dmte$ratio_obs_mod_mte  <-  nice_to_dmte$gpp_obs / nice_to_dmte$gpp_mte
  #       # nice_to_dmte$ratio_obs_mod_mte[ which(is.infinite(nice_to_dmte$ratio_obs_mod_mte)) ] <- NA        

  #       ## save to file
  #       save( nice_to_mte, file=filn )

  #     }

  #     ## add row to aggregated data
  #     mysitename <- data.frame( mysitename=rep( sitename, nrow(nice_to_mte) ) )
  #     nice_to_mte_agg <- rbind( nice_to_mte_agg, cbind( mysitename, dplyr::select( nice_to_mte, bias_mte, is_drought_byvar, gpp_mte, gpp_obs ) ) )

  #     idxs_drought <- which( nice_to_mte$is_drought_byvar )

  #   } else {

  #     missing_mte <- TRUE

  #   }

  # }

  # ##------------------------------------------------
  # ## Get MODIS-GPP for this site
  # ##------------------------------------------------
  # if (avl_data_modis){

  #   filn <- paste( myhome, "data/modis_", sitename, ".Rdata", sep="" )

  #   if ( file.exists(filn) && !overwrite_modis ){

  #     ## load 'nice_to_modis' from file
  #     load( filn )
  #     avl_modisgpp <- TRUE

  #   } else {

  #     ## prepare 'nice_to_modis'
  #     modis <- try( read.csv( paste( myhome, "data/modis_gpp_fluxnet_cutouts_tseries/", sitename, "/gpp_8d_modissubset_", sitename, ".csv", sep="" ), as.is=TRUE ))
  #     if (class(modis)!="try-error"){
  #       avl_modisgpp <- TRUE
  #       modis <- dplyr::rename( modis, gpp_modis=data )
  #       modis$gpp_modis <- modis$gpp_modis / 8

  #       ## Make 'nice' dataframe conform with 'modis'. Take mean of all variables across days for which year_dec is within four +/- 4 days of the modis date
  #       nice_to_modis <- c()
  #       for (idx in 1:nrow(modis)){
  #         sub <- dplyr::filter( nice, year_dec>=(modis$year_dec[idx] - 4/365 ) & year_dec<=(modis$year_dec[idx] + 4/365) )
  #         addline <- unlist( unname( apply( sub, 2, FUN=mean, na.rm=TRUE ) ) )
  #         nice_to_modis <- rbind( nice_to_modis, addline )
  #       }

  #       nice_to_modis <- as.data.frame(nice_to_modis)
  #       colnames( nice_to_modis ) <- names( sub )
  #       rownames( nice_to_modis ) <- NULL

  #       nice_to_modis$year_dec  <- modis$year_dec
  #       nice_to_modis$gpp_modis <- modis$gpp_modis

  #       nice_to_modis$is_drought_byvar <- with( nice_to_modis, ifelse( is_drought_byvar<0.5, FALSE, TRUE ) )

  #       # plot(  nice_to_modis$year_dec, nice_to_modis$gpp_obs, type='l' )
  #       # lines( nice_to_modis$year_dec, nice_to_modis$gpp_modis, col='red' )

  #       # plot( nice_to_modis$year_dec, (nice_to_modis$gpp_modis - nice_to_modis$gpp_obs) / nice_to_modis$gpp_obs, type='l' )

  #       ## get bias measure: log( mod / obs )
  #       nice_to_modis$bias_modis <- nice_to_modis$gpp_modis / nice_to_modis$gpp_obs
  #       nice_to_modis$bias_modis[ which(is.infinite(nice_to_modis$bias_modis)) ] <- NA
  #       nice_to_modis$ratio_obs_mod_modis  <-  nice_to_modis$gpp_obs / nice_to_modis$gpp_modis
  #       nice_to_modis$ratio_obs_mod_modis[ which(is.infinite(nice_to_modis$ratio_obs_mod_modis)) ] <- NA

  #       ## save to file
  #       save( nice_to_modis, file=filn )

  #     } else {
  #       avl_modisgpp <- FALSE
  #     }

  #   }

  #   ## add row to aggregated data
  #   if (avl_modisgpp){
  #     idxs_drought <- which( nice_to_modis$is_drought_byvar )

  #     mysitename <- data.frame( mysitename=rep( sitename, nrow(nice_to_modis) ) )
  #     nice_to_modis_agg <- rbind( nice_to_modis_agg, cbind( mysitename, dplyr::select( nice_to_modis, bias_modis, is_drought_byvar, gpp_modis, gpp_obs ) ) )
  #   }

  # }

}

print("... done.")

if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  save( nice_agg,  file=paste("data/nice_agg_",  nam_target, char_fapar, ".Rdata", sep="") )
  save( nice_resh, file=paste("data/nice_resh_", nam_target, char_fapar, ".Rdata", sep="") )

  if (avl_data_mte)   save( nice_to_mte_agg,   file=paste("data/nice_mte_agg_",   nam_target, char_fapar, ".Rdata", sep="") )
  if (avl_data_modis) save( nice_to_modis_agg, file=paste("data/nice_modis_agg_", nam_target, char_fapar, ".Rdata", sep="") )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

