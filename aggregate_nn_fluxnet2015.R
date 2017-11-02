.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

library( dplyr )
library( cgwtools )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/utils_sofun/analysis_sofun/remove_outliers.R", sep="" ) )

# IMPORTANT: USE SOILMOISTURE FROM S13 FOR NN-TRAINING
load( paste( myhome, "data/fluxnet_sofun/modobs_fluxnet2015_s11_s12_s13_with_SWC_v3.Rdata", sep="" ) ) # "new data" with s13

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read.csv( "successcodes.csv", as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
# do.sites   = "AU-Dry"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
overwrite_nice = FALSE
overwrite_modis = FALSE
overwrite_mte = FALSE
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

print( paste( "Aggregating and complementing data for ", length(do.sites)," NN FLUXNET2015 sites ...") )

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

  nicefiln <- paste0("data/nice_nn_", sitename, ".Rdata" )

  if (file.exists(nicefiln)&&!overwrite_nice){

    load( nicefiln )

  } else {

    infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_", nam_target, char_wgt, char_fapar, ".Rdata", sep="" ) 
    if (verbose) print( paste( "opening file", infil ) )

    load( infil ) ## gets list 'nn_fluxnet'
    nice             <- as.data.frame( nn_fluxnet[[ sitename ]]$nice )            
    varnams_swc      <- nn_fluxnet[[ sitename ]]$varnams_swc    
    varnams_swc_obs  <- nn_fluxnet[[ sitename ]]$varnams_swc_obs

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
    nice <- nice %>%  mutate( bias_pmodel = gpp_pmodel / gpp_obs, 
                              ratio_obs_mod = gpp_obs / gpp_pmodel, 
                              alpha = aet_pmodel / pet_pmodel ) %>% 
                      mutate( bias_pmodel = ifelse( is.infinite(bias_pmodel), NA, bias_pmodel ), 
                              ratio_obs_mod = ifelse( is.infinite(ratio_obs_mod), NA, ratio_obs_mod ),
                              alpha = ifelse( is.infinite(alpha), NA, alpha )
                              )

    ## add row to aggregated data
    nice <- nice %>% mutate( mysitename=sitename )

    ## fLUE estimate based on current soil moisture and average AET/PET
    meanalpha <- mean( nice$aet_pmodel / nice$pet_pmodel, na.rm=TRUE )
    nice <- nice %>%  mutate( dry = ifelse(alpha<0.95, TRUE, FALSE) )

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
    ## save to file
    ##------------------------------------------------
    save( nice, file=nicefiln )

  }    

  ##------------------------------------------------
  ## record for aggregated data
  ##------------------------------------------------
  usecols <- c(
                "mysitename",
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
                "is_drought_byvar", 
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
                "dry",
                paste( "var_nn_act_", varnams_swc_full, sep="" ),
                paste( "var_nn_pot_", varnams_swc_full, sep="" ),
                paste( "var_nn_vpd_", varnams_swc_full, sep="" ),
                paste( "moist_", varnams_swc_full, sep="" )
                )

  nice_agg <- rbind( nice_agg, select( nice, one_of( usecols ) ) )

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

      filn <- paste( "data/mte_nn_", sitename, ".Rdata", sep="" )

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
                          mutate( inmtebin = cut( as.numeric(date), breaks = c( mte$date_start ), right = FALSE ) ) 

        ## summarise by bin taking means
        nice_to_mte <- nice %>% group_by( inmtebin ) %>% summarise_all( mean, na.rm=TRUE )
        tmp         <- nice %>% select( inmtebin, doy, year ) %>% group_by( inmtebin ) %>% summarise( doy_start=min( doy ), doy_end=max( doy ), year_start=min( year ), year_end=max( year ) )
        nice_to_mte <- cbind( nice_to_mte, select( tmp, doy_start, doy_end, year_start, year_end ) )

        ## merge dataframes (averaged nice and mte)
        nice_to_mte <- nice_to_mte %>% select( doy_start, year_start, one_of( usecols ), inmtebin, -mysitename ) %>% left_join( mte, by=c("doy_start", "year_start") )
        
        ## get additional variables
        nice_to_mte <- nice_to_mte %>%  mutate( bias_mte = gpp_mte / gpp_obs )          %>% mutate( bias_mte=ifelse( is.infinite( bias_mte ), NA, bias_mte ) ) %>% 
                                        mutate( ratio_obs_mod_mte = gpp_obs / gpp_mte ) %>% mutate( ratio_obs_mod_mte=ifelse( is.infinite( bias_mte ), NA, ratio_obs_mod_mte ) ) %>%  
                                        mutate( bias_rf = gpp_rf / gpp_obs )            %>% mutate( bias_rf=ifelse( is.infinite( bias_rf ), NA, bias_rf ) )  %>%  
                                        mutate( ratio_obs_mod_rf = gpp_obs / gpp_rf )   %>% mutate( ratio_obs_mod_rf=ifelse( is.infinite( bias_rf ), NA, ratio_obs_mod_rf ) ) %>% 
                                        mutate( is_drought_byvar = ifelse( is_drought_byvar>=0.5, TRUE, FALSE ) )

        ## save to file
        save( nice_to_mte, file=filn )

      }

      ## add row to aggregated data
      mte_agg <- rbind( mte_agg, select( nice_to_mte, one_of( c( usecols, "bias_mte", "ratio_obs_mod_mte", "bias_rf", "ratio_obs_mod_rf", "doy_start", "doy_end", "year_start", "year_end" ) ) ) )

    } else {

      missing_mte <- TRUE

    }

  }


  ##------------------------------------------------
  ## Get MODIS-GPP for this site
  ##------------------------------------------------
  if (avl_data_modis){

    filn <- paste( "data/modis_nn_", sitename, ".Rdata", sep="" )

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
        modis <- modis %>%  rename( gpp_modis = data ) %>% 
                            mutate( gpp_modis = gpp_modis / 8.0, doy_start = doy, doy_end = lead( doy ) - 1, year_start = year, 
                                    date_start = as.POSIXct( as.Date( date ) ), date_end = as.POSIXct( as.Date( lead( date ) ) - 1 ),
                                    mysitename = sitename
                                    )

        ## group nice by 8d bins from MODIS data
        nice <- nice %>%  mutate( date = as.POSIXct( as.Date( paste( as.character(year), "-01-01", sep="" ) ) + doy - 1 ) ) %>% 
                          mutate( inmodisbin = cut( as.numeric(date), breaks = c( modis$date_start ), right = FALSE ) )

        ## summarise by bin taking means
        nice_to_modis <- nice %>% group_by( inmodisbin ) %>% summarise_all( mean, na.rm=TRUE )
        tmp           <- nice %>% select( inmodisbin, doy, year ) %>% group_by( inmodisbin ) %>% summarise( doy_start=min( doy ), doy_end=max( doy ), year_start=min( year ), year_end=max( year ) )
        nice_to_modis <- cbind( nice_to_modis, select( tmp, doy_start, doy_end, year_start, year_end ) )

        ## merge dataframes (averaged nice and modis)
        nice_to_modis <- nice_to_modis %>% select( doy_start, year_start, one_of( usecols ), inmodisbin, -mysitename ) %>% left_join( modis, by=c("doy_start", "year_start") )
        
        ## get additional variables
        nice_to_modis <- nice_to_modis %>%  mutate( bias_modis = gpp_modis / gpp_obs )          %>% mutate( bias_modis=ifelse( is.infinite( bias_modis ), NA, bias_modis ) ) %>% 
                                            mutate( ratio_obs_mod_modis = gpp_obs / gpp_modis ) %>% mutate( ratio_obs_mod_modis=ifelse( is.infinite( bias_modis ), NA, ratio_obs_mod_modis ) ) %>% 
                                            mutate( is_drought_byvar = ifelse( is_drought_byvar>=0.5, TRUE, FALSE ) )

        ## save to file
        save( nice_to_modis, file=filn )

      } else {

        avl_modisgpp <- FALSE

      }

    }

    ## add row to aggregated data
    if (avl_modisgpp){
      modis_agg <- rbind( modis_agg, select( nice_to_modis, one_of( c( usecols, "bias_modis", "ratio_obs_mod_modis", "doy_start", "doy_end", "year_start", "year_end" )) ) )
    }

  }

}

print("... done.")
print("Aggregation completed.")
print( paste( "site data with P-model outputs has been saved in files  data/nice_nn_<sitename>.Rdata"))
print( paste( "site data with FLUXCOM MTE data has been saved in files data/mte_nn_<sitename>.Rdata"))
print( paste( "site data with MODIS GPP data has been saved in files   data/modis_nn_<sitename>.Rdata"))


if ( length( dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename ) == length( do.sites ) ){
  ##------------------------------------------------
  ## save collected data
  ##------------------------------------------------
  filn <- paste0("data/nice_nn_agg_",  nam_target, char_fapar, ".Rdata")
  print( paste( "saving dataframe nice_agg in file", filn) )
  save( nice_agg,  file=filn )

  filn <- paste0("data/nice_nn_mte_agg_",   nam_target, char_fapar, ".Rdata")
  print( paste( "saving dataframe mte_agg in file", filn) )
  save( mte_agg, file=filn )

  filn <- paste0("data/nice_nn_modis_agg_", nam_target, char_fapar, ".Rdata")
  print( paste( "saving dataframe modis_agg in file", filn) )
  save( modis_agg, file=filn )

  # save( nice_resh, file="data/nice_resh_lue_obs_evi.Rdata" )
  # save( overview, file="data/overview_data_fluxnet2015_L3.Rdata" )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

