library(dplyr)
library(R.matlab)
library(readr)
library(lubridate)
library(tibble)

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

source( paste( myhome, "sofun/utils_sofun/analysis_sofun/remove_outliers.R", sep="" ) )
source( "../utilities/init_dates_dataframe.R" )
source("calc_flue_est_alpha.R")
source("stress_quad_1sided.R")
source("stress_exp.R")

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read_csv( "successcodes.csv" )

# ## Exclude sites for which no fapar data is available and hence no model results
# df_error_fapar <- read_csv( paste0( myhome, "sofun/input_fluxnet2015_sofun/error_missing_forcingdata_MODIS_FPAR_MCD15A3H_fluxnet2015.csv" ) ) 
# successcodes <- successcodes %>% left_join( df_error_fapar, by="mysitename" ) %>% rename( error_fapar = error ) %>% filter( error_fapar == 0 )

do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Manual settings ----------------
# do.sites   = "FR-Pue" # uncomment to run for single site
simsuite = "fluxnet2015"
outputset = "s16"
nam_target = "lue_obs_evi"
use_weights= FALSE    
use_fapar  = FALSE
package    = "nnet"
overwrite_nice = TRUE
overwrite_8d    = TRUE
verbose         = FALSE
##---------------------------------

##------------------------------------------------
## Get FLUXNET 2015 data and SOFUN outputs from site-scale simulations
## The file loaded here is created by 'get_modobs.R'
##------------------------------------------------
datafilnam_flat <- paste0( "data/df_modobs_fluxnet2015_", paste( outputset, collapse="_" ), "_with_SWC_v5.Rdata" )
load( datafilnam_flat )  # loads 'df_fluxnet'

##------------------------------------------------
## Load data for NN predictions
##------------------------------------------------
## load the NN model object
# load( "data/nn_caret_GRA.Rdata" )
# load( "data/rf.Rdata" )
load( "data/ai_fluxnet2015.Rdata" )  # loads 'df_ai' (aridity index)
load( "data/wtd_fluxnet2015.Rdata" )  # loads 'df_wtd' (water table depth)
siteinfo <- read_csv( "data/siteinfo_fluxnet2015_sofun.csv" ) # needed for vegetation type

##------------------------------------------------
## Get MTE-GPP for all sites
##------------------------------------------------
filn <- paste0( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_8Days_4Beni.csv" )
if ( file.exists( filn ) ){
  mte_dl <- read_csv( paste0( myhome, "data/gpp_mte_rf_fluxnet_tramontana/GPP_Daily_4Beni.csv" ) )
  mte_8d <- read_csv( filn ) %>% rename(  mysitename=`Site-code`, 
                                          gpp_mte=MTE, 
                                          gpp_mte_m=MTE_M, 
                                          gpp_mte_viterbo=MTE_Viterbo, 
                                          gpp_rf=RF, 
                                          gpp_rf_fromdaily=RF_from_daily, 
                                          doy_start=StartDoY, 
                                          doy_end=EndDoY, 
                                          year_start=StartYear, 
                                          year_end=EndYear 
                                          ) %>%
                                  mutate( mysitename=as.character(mysitename),
                                          year_start=as.integer(year_start),
                                          year_end=as.integer(year_end),
                                          doy_start=as.integer(doy_start),
                                          doy_end=as.integer(doy_end)
                                          )
}

## replace with NA
for (ivar in names(mte_8d)){ mte_8d[[ ivar ]][ which( mte_8d[[ ivar ]]==-9999 ) ] <- NA }


##------------------------------------------------
## Get BESS-GPP for all sites
##------------------------------------------------
if ( file.exists( paste0( myhome, "data/gpp_bess/sitescale_fluxnet/BESSv1.GPP.Daily.mat" ) ) ){

  df_bess_gpp_v1 <- readMat( paste0( myhome, "data/gpp_bess/sitescale_fluxnet/BESSv1.GPP.Daily.mat" ) )$data %>% as.data.frame()
  df_bess_gpp_v2 <- readMat( paste0( myhome, "data/gpp_bess/sitescale_fluxnet/BESSv2.GPP.Daily.mat" ) )$data %>% as.data.frame()

  meta <- read_csv( paste0( myhome, "data/gpp_bess/sitescale_fluxnet/FLUXNET2015v1-3.csv") )

  colnames(df_bess_gpp_v1) <- meta$Name
  colnames(df_bess_gpp_v2) <- meta$Name

  df_tmp <- init_dates_dataframe( 2001, 2014, noleap=TRUE )

  df_bess_gpp_v1 <- df_tmp %>% bind_cols( df_bess_gpp_v1 )
  df_bess_gpp_v2 <- df_tmp %>% bind_cols( df_bess_gpp_v2 )

  avl_data_bess <- TRUE

} else {

  avl_data_bess <- FALSE

}

##------------------------------------------------
## Get VPM-GPP for all sites
##------------------------------------------------
filn <- paste0( myhome, "data/gpp_vpm/sitescale_fluxnet/gpp_vpm_combined_data_for_figure_allsites_v20.csv" )
if ( file.exists( filn ) ){

  df_gpp_vpm <- read_csv( filn ) %>% 
                mutate( year=as.integer( substr( as.character(date), start=1, stop=4 )), 
                        doy=as.integer( substr( as.character(date), start=5, stop=7 ) ) ) %>%
                rename( mysitename=sid ) %>%
                select( -n, -biome, -date ) %>%
                mutate( date = ymd( paste0( as.character(year), "-01-01" ) ) + days(doy - 1) )

}

print( paste( "Aggregating and complementing data for", length(do.sites), "NN FLUXNET2015 sites ...") )

##------------------------------------------------
## Initialise aggregated data
##------------------------------------------------
## fvar and soilm data to be complemented with cluster info
nice_agg  <- data.frame()
nice_8d_agg <- data.frame()

for (sitename in do.sites){

  print( paste( sitename, "..."))

  missing_mte <- FALSE

  nicefiln <- paste0("data/nice_nn/nice_nn_", sitename, ".Rdata" )

  if (file.exists(nicefiln)&&!overwrite_nice){

    load( nicefiln )

  } else {

    ##------------------------------------------------
    ## load site data and "detatch"
    ##------------------------------------------------
    nice <- filter( df_fluxnet, mysitename==sitename ) %>% 
            select( date, temp, ppfd, vpd, prec, fpar, gpp_obs=GPP_NT_VUT_REF, starts_with("SWC_F_MDS"), soilm_splash, gpp_pmodel=gpp, aet_pmodel=aet, pet_pmodel=pet )

    ##------------------------------------------------
    ## Get NN-derived data (flue, is_drought_byvar) from separate file
    ##------------------------------------------------
    infil <- paste( myhome, "data/nn_fluxnet/fvar/nn_fluxnet2015_", sitename, "_lue_obs_evi.Rdata", sep="" ) 
    load( infil ) ## gets list 'nn_fluxnet'
    nn_nice <-  as_tibble( nn_fluxnet[[ sitename ]]$nice ) %>%
                mutate( date = ymd( paste0( as.character(year), "-01-01" ) ) + days(doy - 1) ) %>%
                select( -gpp_pmodel, -temp, -ppfd, -vpd, -prec, -evi, -fpar, -lue_obs_evi, -lue_obs_fpar, -aet_pmodel, -pet_pmodel, -soilm_splash150, -soilm_splash220, -soilm_swbm, -soilm_etobs, -soilm_etobs_ob, -gpp_obs )

    ##------------------------------------------------
    ## nice and nn_nice have common dates. merge them.
    ## Make sure to keep variables that are now already in 'nice'.
    ##------------------------------------------------
    nice <- nn_nice %>% left_join( nice, by="date")
    
    ##------------------------------------------------
    ## additional variables
    ##------------------------------------------------
    nice <- nice %>%  mutate( bias_pmodel = gpp_pmodel / gpp_obs,
                              bias_pmodel_diff = gpp_pmodel - gpp_obs,
                              ratio_obs_mod_pmodel = gpp_obs / gpp_pmodel, 
                              alpha = aet_pmodel / pet_pmodel ) %>% 

                      mutate( bias_pmodel = ifelse( is.infinite(bias_pmodel), NA, bias_pmodel ), 
                              ratio_obs_mod_pmodel = ifelse( is.infinite(ratio_obs_mod_pmodel), NA, ratio_obs_mod_pmodel ),
                              alpha = ifelse( is.infinite(alpha), NA, alpha )
                              ) %>%

                      ## add site name as column
                      mutate( mysitename=sitename ) %>%

                      ## date as ymd from the lubridate package
                      mutate( date = ymd( paste0( as.character(year), "-01-01" ) ) + days( doy - 1 ) ) %>%

                      ## dry
                      mutate( dry = ifelse(alpha<0.95, TRUE, FALSE) ) %>%

                      ## get LUE and remove outliers
                      mutate( lue_obs_fpar = remove_outliers( gpp_obs / ( ppfd * fpar ), coef=3.0 ) ) %>%

                      ## get mean observational soil moisture across different depths (if available)
                      mutate( soilm_obs_mean = apply( select( ., starts_with("SWC_F_MDS") ), 1, FUN = mean, na.rm = TRUE ) ) %>%
                      mutate( soilm_obs_mean = ifelse( is.nan(soilm_obs_mean), NA, soilm_obs_mean ) )



    ## fLUE estimate based on current soil moisture and average AET/PET
    meanalphaval <- mean( nice$aet_pmodel / nice$pet_pmodel, na.rm=TRUE )

    ## add vegetation type to dataframe
    nice <- nice %>% left_join( select( siteinfo, mysitename, classid ), by="mysitename")

    ##------------------------------------------------
    ## Calculate soil moiture stress based on two alternative functions
    ##------------------------------------------------
    ## soil moisture stress factor (derived from two different fitting methods, see knit...Rmd)
    myclassid <- unique( nice$classid )
    mymeanalpha <- meanalphaval
    nice <- nice %>%  mutate( meanalpha=meanalphaval ) %>%
                      mutate( flue_est_I  = calc_flue_est_alpha( soilm_mean, mymeanalpha, apar=0.2617121, bpar=0.5668587, cpar=0.125, dpar=0.75 ), ## method for s1a
                              flue_est_II = stress_quad_1sided_alpha( soilm_mean, mymeanalpha, x0 = 0.9, apar = 0.0606651, bpar = 0.5090085 ),  ## when fitting to directly to ratio_obs_mod_pmodel, subset of notoriously sensitive sites
                              flue_est_III = stress_quad_1sided_alpha( soilm_mean, mymeanalpha, x0 = 0.9, apar = 0.05151085, bpar = 0.1920844 )  ## when fitting to directly to ratio_obs_mod_pmodel, subset of notoriously sensitive sites
                              )

    nice$flue_est_IV <- rep( NA, nrow( nice ))
    nice$flue_est_V  <- rep( NA, nrow( nice ))
    for (idx in seq(nrow(nice))){
      nice$flue_est_IV[idx] <- stress_quad_1sided_alpha_grasstree( nice$soilm_mean[idx], mymeanalpha, x0 = 0.9, apar = c( 0.1785247, 0.1007749 ), bpar = c(0.4501654, 0.0063069), classid=myclassid )  ## when fitting to ratio_obs_mod_pmodel, method for s1b
      nice$flue_est_V[idx]  <- stress_exp_alpha_grasstree( nice$soilm_mean[idx], mymeanalpha, c(0.2067347, 0.0056745), c( 0.5543440, 0.4116154), 6.170356, -3.797436, classid=myclassid ) ## when fitting to directly to ratio_obs_mod_pmodel, subset of notoriously sensitive sites)
    }


    # ##------------------------------------------------
    # ## Add NN-based fLUE estimate
    # ##------------------------------------------------
    # ## add mean aridity index (P/PET) to dataframe
    # nice <- nice %>% left_join( select( df_ai, mysitename, ai ), by="mysitename" )

    # ## add water table depth to dataframe
    # nice <- nice %>% left_join( df_wtd, by="mysitename" )

    # ## predict the values
    # tmp <- select( nice, mysitename, date, soilm_splash, classid, ai, wtd, fpar) %>%
    #        mutate( GRA = ifelse( classid=="GRA", TRUE, FALSE ),
    #                SAV = ifelse( classid=="SAV", TRUE, FALSE ),
    #                ENF = ifelse( classid=="ENF", TRUE, FALSE ),
    #                WET = ifelse( classid=="WET", TRUE, FALSE ),
    #                WSA = ifelse( classid=="WSA", TRUE, FALSE ),
    #                EBF = ifelse( classid=="EBF", TRUE, FALSE ),
    #                DBF = ifelse( classid=="DBF", TRUE, FALSE ),
    #                CSH = ifelse( classid=="CSH", TRUE, FALSE )
    #               ) %>%
    #        mutate(  classid = as.factor(classid), 
    #                 GRA = as.factor(GRA),
    #                 SAV = as.factor(SAV),
    #                 ENF = as.factor(ENF),
    #                 WET = as.factor(WET),
    #                 WSA = as.factor(WSA),
    #                 EBF = as.factor(EBF),
    #                 DBF = as.factor(DBF),
    #                 CSH = as.factor(CSH)
    #                 )
    # na_idxs <- apply( tmp, 2, FUN=function (x) which(is.na(x)) ) %>% unlist() %>% unique()
    # tmp <- tmp[-na_idxs,]
    # print(  apply( tmp, 2, FUN=function (x) sum(is.na(x)) ) )
    
    # ## predict
    # vals <- as.vector( predict( rf, tmp ) )
    # tmp <- tmp %>% mutate( flue_est_nn = vals )
    # nice <- nice %>% left_join( select( tmp, mysitename, date, flue_est_nn ), by = c("mysitename", "date") )

    # ## test plot
    # with(nice, plot(date, fvar, type="l", main=sitename, ylim=c(0,1)))
    # with(nice, lines(date, flue_est_nn, col="red"))
    
    ##------------------------------------------------
    ## Add BESS-GPP for this site to daily dataframe 'nice'
    ##------------------------------------------------
    if (avl_data_bess){

      if (is.element( sitename, colnames(df_bess_gpp_v1))){

          print( "getting BESS data ...")

          ## get bess v1 data for this site
          bess <- df_bess_gpp_v1[ , is.element( names(df_bess_gpp_v1), c("date", sitename) )] %>%
                    setNames( c("date", "gpp_bess_v1")) 

          ## add bess v2 data for this site                                      
          bess <- df_bess_gpp_v2[ , is.element( names(df_bess_gpp_v2), c("date", sitename) )] %>%
                    setNames( c("date", "gpp_bess_v2")) %>%
                    right_join( bess, by="date" )
        
          ## merge into nice dataframe
          nice <- nice %>%  left_join( bess, by = "date" ) %>% 
                            mutate( bias_bess_v1 = gpp_bess_v1 / gpp_obs, bias_bess_v1_diff = gpp_bess_v1 - gpp_obs ) %>% mutate( bias_bess_v1 = ifelse( is.infinite( bias_bess_v1 ), NA, bias_bess_v1 ) ) %>% 
                            mutate( ratio_obs_mod_bess_v1 = gpp_obs / gpp_bess_v1 ) %>% mutate( ratio_obs_mod_bess_v1=ifelse( is.infinite( bias_bess_v1 ), NA, ratio_obs_mod_bess_v1 ) )  %>%
                            mutate( bias_bess_v2 = gpp_bess_v2 / gpp_obs, bias_bess_v2_diff = gpp_bess_v2 - gpp_obs ) %>% mutate( bias_bess_v2 = ifelse( is.infinite( bias_bess_v2 ), NA, bias_bess_v2 ) ) %>% 
                            mutate( ratio_obs_mod_bess_v2 = gpp_obs / gpp_bess_v2 ) %>% mutate( ratio_obs_mod_bess_v2=ifelse( is.infinite( bias_bess_v2 ), NA, ratio_obs_mod_bess_v2 ) )

        }

    }

    ##------------------------------------------------
    ## save to file
    ##------------------------------------------------
    if (verbose) print( paste( "saving to file", nicefiln ) )
    save( nice, file=nicefiln )

  }


  ##------------------------------------------------
  ## record for aggregated data
  ##------------------------------------------------
  usecols <- c(
                "mysitename",
                "date",
                "gpp_obs",
                "var_nn_pot", 
                "var_nn_act",
                "var_nn_vpd",
                "gpp_nn_act", 
                "gpp_nn_pot",
                "fvar",
                "fvar_smooth",
                "flue_est_I",
                "flue_est_II",
                "flue_est_III",
                "flue_est_IV",
                "flue_est_V",
                "is_drought_byvar", 
                "gpp_pmodel",
                "aet_pmodel",
                "pet_pmodel",
                "alpha",
                "ppfd",
                "vpd",
                "evi", 
                "fpar", 
                "soilm_splash",
                "soilm_swbm",
                "soilm_mean",
                "soilm_obs_mean",
                "bias_pmodel", 
                "bias_pmodel_diff", 
                "ratio_obs_mod_pmodel", 
                "lue_obs_evi", 
                "lue_obs_fpar",
                "dry",
                "gpp_bess_v1", "gpp_bess_v2",
                "bias_bess_v1", "bias_bess_v2",
                "bias_bess_v1_diff", "bias_bess_v2_diff",
                "ratio_obs_mod_bess_v1", "ratio_obs_mod_bess_v2"
                )

  nice_agg <- bind_rows( nice_agg, select( nice, one_of( usecols ) ) )


  ##------------------------------------------------
  ## Aggregate 'nice' to 8-day periods from the MODIS data and get MODIS-GPP for this site
  ##------------------------------------------------
  filn_8d <- paste0( "data/nice_8d_nn/nice_8d_nn", sitename, ".Rdata" )
  if (!dir.exists("data/nice_8d_nn")) system("mkdir -p data/nice_8d_nn")

  if ( file.exists(filn_8d) && !overwrite_8d ){

    ## load 'nice_8d' from file
    load( filn_8d )

  } else {

    print("getting MODIS data")

    ## prepare 'nice_8d'
    modis <- try( read_csv( paste0( myhome, "data/gpp_MODIS_GPP_MOD17A2H_fluxnet2015_gee_subset/gpp_MODIS_GPP_MOD17A2H_", sitename, "_gee_subset.csv" ) ) )

    if (class(modis)[1]!="try-error"){

      avl_modisgpp <- TRUE

      ## make modis a bit nicer
      modis <- modis %>%  rename( gpp_modis = modisvar ) %>% 
                          mutate( gpp_modis = ifelse( !is.na(gpp_modis), gpp_modis * 1e3 / 8.0, NA ),
                                  date_end = ymd( lead( date_start ) ) - days(1) ) %>%
                          select( date_start, date_end, gpp_modis )

      ## group nice by 8d bins from MODIS data
      nice <- nice %>%  mutate( in8dbin = cut( date, breaks = modis$date_start, right = FALSE ) )

      ## summarise by bin taking means
      nice_8d <- nice %>% group_by( in8dbin ) %>% summarise_all( mean, na.rm=TRUE )
      tmp     <- nice %>% select( in8dbin, date ) %>% group_by( in8dbin ) %>% summarise( date_start=min( date ), date_end=max( date ) )
      nice_8d <- bind_cols( nice_8d, select( tmp, date_start, date_end ) )

      # print( modis$date_start %in% nice_8d$date_start )  
      
      ## merge modis data into 8-day dataframe (averaged nice and modis)
      nice_8d <- nice_8d %>% select( date_start, one_of( usecols ) ) %>% left_join( select( modis, -date_end ), by=c( "date_start" ) )
      
      ## get additional variables
      nice_8d <- nice_8d %>% mutate( bias_modis = gpp_modis / gpp_obs, bias_modis_diff = gpp_modis - gpp_obs ) %>% mutate( bias_modis=ifelse( is.infinite( bias_modis ), NA, bias_modis ) ) %>% 
                             mutate( ratio_obs_mod_modis = gpp_obs / gpp_modis ) %>% mutate( ratio_obs_mod_modis=ifelse( is.infinite( ratio_obs_mod_modis ), NA, ratio_obs_mod_modis ) ) %>%
                             mutate( mysitename = sitename )


      ##------------------------------------------------
      ## Get MTE-GPP for this site
      ##------------------------------------------------
      if ( is.element( sitename, mte_8d$mysitename ) ){

        print( "getting MTE data ...")

        ## prepare dataframe 'nice_to_mte'
        # dmte <- filter( mte_dl, Site.code==sitename )
        avl_mtegpp <- TRUE

        ## make mte a bit nicer
        mte <-  filter( mte_8d, mysitename==sitename ) %>%
                mutate( doy_end = doy_end - 1 ) %>%  # assuming that the doy_end is not counted towards bins aggregate
                mutate( date_start = ymd( paste0( as.character(year_start), "-01-01" ) ) + days( doy_start - 1 ),
                        date_end   = ymd( paste0( as.character(year_end  ), "-01-01" ) ) + days( doy_end   - 1 ) ) %>%
                select( -ID, -year_start, -doy_start, -year_end, -doy_end )

        # ## check if the bins in the MTE data are identical to the ones in the MODIS data
        # print( mte$date_start %in% nice_8d$date_start )  

        ## merge to MTE data to nice_8d
        nice_8d <- nice_8d %>% left_join( select( mte, -date_end ), by = c("date_start", "mysitename") ) %>%

          ## get additional variables
          mutate( bias_mte = gpp_mte / gpp_obs, bias_mte_diff = gpp_mte - gpp_obs ) %>% mutate( bias_mte=ifelse( is.infinite( bias_mte ), NA, bias_mte ) ) %>% 
          mutate( ratio_obs_mod_mte = gpp_obs / gpp_mte )                           %>% mutate( ratio_obs_mod_mte=ifelse( is.infinite( ratio_obs_mod_mte ), NA, ratio_obs_mod_mte ) ) %>%  
          mutate( bias_rf = gpp_rf / gpp_obs, bias_rf_diff = gpp_rf - gpp_obs )     %>% mutate( bias_rf=ifelse( is.infinite( bias_rf ), NA, bias_rf ) )  %>%  
          mutate( ratio_obs_mod_rf = gpp_obs / gpp_rf )                             %>% mutate( ratio_obs_mod_rf=ifelse( is.infinite( ratio_obs_mod_rf ), NA, ratio_obs_mod_rf ) )

      } else {

        avl_mtegpp <- FALSE

      }


      ##------------------------------------------------
      ## Get VPP-GPP for this site
      ##------------------------------------------------    
      vpm <- df_gpp_vpm %>% filter( mysitename==sitename )
      if (nrow(vpm)>0){

        print( "getting VPM data ...")

        avl_vpmgpp <- TRUE

        ## make vpm a bit nicer
        vpm <- vpm %>%  select( -GPP ) %>% rename( gpp_vpm = VPM ) %>% 
                        mutate( date_start = date, date_end = lead( date ) - days(1) )

        # ## check if the bins in the MTE data are identical to the ones in the MODIS data
        # print( is.element( vpm$date_start, nice_8d$date_start ) )

        ## merge to MTE data to nice_8d
        nice_8d <- nice_8d %>% left_join( select( vpm, -date, -date_end, -year, -doy ), by = c("date_start", "mysitename") ) %>%

          ## get additional variables
          mutate( bias_vpm = gpp_vpm / gpp_obs, bias_vpm_diff = gpp_vpm - gpp_obs ) %>% mutate( bias_vpm=ifelse( is.infinite( bias_vpm ), NA, bias_vpm ) ) %>% 
          mutate( ratio_obs_mod_vpm = gpp_obs / gpp_vpm ) %>% mutate( ratio_obs_mod_vpm=ifelse( is.infinite( ratio_obs_mod_vpm ), NA, ratio_obs_mod_vpm ) )
      
      } else {

        avl_vpmgpp <- FALSE

      }


    } else {

      avl_modisgpp <- FALSE

    }

  }

  ## save to file
  if (verbose) print( paste( "saving to file", filn_8d ) )
  save( nice_8d, file=filn_8d )

  ## add row to aggregated data
  nice_8d_agg <- bind_rows( nice_8d_agg, nice_8d )

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
  filn <- paste0("data/nice_nn_agg_lue_obs_evi.Rdata")
  print( paste( "saving dataframe nice_agg in file", filn) )
  save( nice_agg,  file=filn )

  ## save aggregated NN mte
  filn <- paste0("data/nice_nn_8d_agg_lue_obs_evi.Rdata")
  print( paste( "saving dataframe mte_agg in file", filn) )
  save( nice_8d_agg, file=filn )

} else {

  print("WARNING: NO SAVING AT THE END!")

}

