reshape_align_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", bysm=FALSE, use_fapar=FALSE, use_weights=FALSE, overwrite=TRUE, verbose=FALSE ){

  # ## debug-------------------
  # # sitename = "AU-Dry"
  # nam_target="lue_obs_evi"
  # bysm=FALSE
  # use_fapar=FALSE
  # use_weights=FALSE
  # overwrite=TRUE
  # verbose=TRUE
  # #--------------------------

  .libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  require( dplyr )
  require( tidyr )
  #require( cgwtools )

  source( paste( myhome, "sofun/utils_sofun/analysis_sofun/get_consecutive.R", sep="" ) ) 

  ## check and override if necessary
  if ( nam_target=="lue_obs" || nam_target=="lue_obs_evi" || nam_target=="lue_obs_fpar" ){
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

  if (bysm){
    char_bysm <- "_bysm"
  } else {
    char_bysm <- ""
  }

  before <- 30
  after  <- 300

  ## Bins for different variables
  fvarbins  <- seq( from=-20, to=40, by=20 )
  faparbins <- seq( from=-20, to=40, by=20 )
  iwuebins  <- seq( from=-30, to=60, by=30 )

  bincentres_fvar  <- fvarbins[1:(length(fvarbins)-1)]   + (fvarbins[2]-fvarbins[1])/2
  bincentres_fapar <- faparbins[1:(length(faparbins)-1)] + (faparbins[2]-faparbins[1])/2
  bincentres_iwue  <- iwuebins[1:(length(iwuebins)-1)]   + (iwuebins[2]-iwuebins[1])/2

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
                "dry"
                # "flue_est"
                )

  ##------------------------------------------------
  ## load "nice" data for this site
  ##------------------------------------------------
  if (verbose) print("loading nice file ...")
  infil <- paste0( "data/nice_nn/nice_nn_", sitename, ".Rdata" )
  load( infil )
  df <- nice; rm("nice")

  ##------------------------------------------------
  ## load drought data (only!) for this site from NN FLUXNET2015 output
  ##------------------------------------------------
  droughtfil <- paste0( "./data/droughts/droughts_", sitename, ".Rdata" )
  if (!file.exists(droughtfil)){
    dir <- paste( myhome, "data/nn_fluxnet/fvar/", sep="" )
    infil <- paste( dir, "nn_fluxnet2015_", sitename, "_", nam_target, char_fapar, ".Rdata", sep="" ) 
    if (verbose) print( paste( "reading file", infil ) )
    load( infil )     ## gets list 'nn_fluxnet'
    droughts <- nn_fluxnet[[ sitename ]]$droughts
    save( droughts, file=droughtfil )    
  } else {
    load( droughtfil )
  }

  names_alg <- c( names(df), "dday")

  ##--------------------------------------------------------
  ## re-arrange data, aligning by beginning of drought
  ##--------------------------------------------------------
  if (nrow(droughts)>1){

    filn <- paste( "./data/df_dday/df_dday_", sitename, ".Rdata", sep="" )

    if (!file.exists(filn)||overwrite){
      if (verbose) print( paste( "aligning df ", sitename, "..." ) )

      df_dday <- data.frame()
      df <- df %>% mutate( mysitename=sitename ) %>% select( one_of(usecols) )
      for ( iinst in 1:nrow(droughts) ){
        after_inst <- min( after, droughts$len[iinst] )
        dday <- seq( from=-before, to=after_inst, by=1 )
        idxs <- dday + droughts$idx_start[iinst]
        drophead <- which( idxs < 1 )
        if (length(drophead)>0){
          idxs <- idxs[ -drophead ]
          dday <- dday[ -drophead ]
        }
        addrows <- df %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst )
        df_dday <- rbind( df_dday, addrows )
      }

      ##--------------------------------------------------------
      ## Bin aligned data and expand from 3D array to dataframe
      ##--------------------------------------------------------
      ## add bin information based on dday to expanded df
      df_dday <- df_dday %>% mutate( infvarbin  = cut( as.numeric(dday), breaks = fvarbins ) )
      df_dday <- df_dday %>% mutate( infaparbin = cut( as.numeric(dday), breaks = faparbins ) )
      df_dday <- df_dday %>% mutate( iniwuebin  = cut( as.numeric(dday), breaks = iwuebins ) )

      ## add row: normalised VPD
      tmp <- df_dday %>% group_by( infvarbin ) %>% 
                         summarise( vpd  = median( vpd , na.rm=TRUE ) ) %>%
                         complete( infvarbin, fill = list( vpd  = NA ) ) %>% 
                         dplyr::select( vpd )
      tmp <- unlist( tmp )[1:(length(fvarbins)-1)]
      df_dday$dvpd = df_dday$vpd / tmp[1]

      ## add row: normalised soil moisture
      tmp <- df_dday %>% group_by( infvarbin ) %>% 
                         summarise( soilm_mean  = median( soilm_mean , na.rm=TRUE ) ) %>%
                         complete( infvarbin, fill = list( soilm_mean  = NA ) ) %>% 
                         dplyr::select( soilm_mean )
      tmp <- unlist( tmp )[1:(length(fvarbins)-1)]
      df_dday$soilm_norm = df_dday$soilm_mean / tmp[1]

      ## add row: normalised fAPAR (EVI)
      tmp <- df_dday %>% group_by( infaparbin ) %>% 
                         summarise( evi  = median( evi , na.rm=TRUE ) ) %>%
                         complete( infaparbin, fill = list( evi  = NA ) ) %>% 
                         dplyr::select( evi )
      tmp <- unlist( tmp )[1:(length(faparbins)-1)]
      df_dday$evi_norm = df_dday$evi / tmp[1]


      ## aggregate by 'dday'
      df_dday_aggbydday <- df_dday %>%  group_by( dday ) %>% 
                                        summarise(
                                                  ## soil moisture
                                                  soilm_med=median( soilm_mean, na.rm=TRUE ), soilm_upp=quantile( soilm_mean, 0.75, na.rm=TRUE ), soilm_low=quantile( soilm_mean, 0.25, na.rm=TRUE ),
                                                  soilm_norm_med=median( soilm_norm, na.rm=TRUE ), soilm_norm_upp=quantile( soilm_norm, 0.75, na.rm=TRUE ), soilm_norm_low=quantile( soilm_norm, 0.25, na.rm=TRUE ),

                                                  ## VPD
                                                  vpd_med=median( vpd, na.rm=TRUE ), vpd_upp=quantile( vpd, 0.75, na.rm=TRUE ), vpd_low=quantile( vpd, 0.25, na.rm=TRUE ),

                                                  ## relative VPD change
                                                  dvpd_med=median( dvpd, na.rm=TRUE ), dvpd_upp=quantile( dvpd, 0.75, na.rm=TRUE ), dvpd_low=quantile( dvpd, 0.25, na.rm=TRUE ),

                                                  ## fLUE
                                                  fvar_med=median( fvar, na.rm=TRUE ), fvar_upp=quantile( fvar, 0.75, na.rm=TRUE ), fvar_low=quantile( fvar, 0.25, na.rm=TRUE ),

                                                  ## EVI and FPAR
                                                  evi_med=median( evi, na.rm=TRUE ), evi_upp=quantile( evi, 0.75, na.rm=TRUE ), evi_low=quantile( evi, 0.25, na.rm=TRUE ),
                                                  fpar_med=median( fpar, na.rm=TRUE ), fpar_upp=quantile( fpar, 0.75, na.rm=TRUE ), fpar_low=quantile( fpar, 0.25, na.rm=TRUE ),

                                                  ## P-model bias
                                                  bias_pmodel_med=median( bias_pmodel, na.rm=TRUE ), bias_pmodel_upp=quantile( bias_pmodel, 0.75, na.rm=TRUE ), bias_pmodel_low=quantile( bias_pmodel, 0.25, na.rm=TRUE )
                                                  ) %>%
                                        mutate( mysitename=sitename )

      ## drop rows dday=NA
      df_dday_aggbydday <- df_dday_aggbydday[ which( !is.na(df_dday_aggbydday$dday)), ]
      save( df_dday_aggbydday, file=paste( "data/df_dday_aggbydday/df_dday_aggbydday_", sitename, ".Rdata", sep="" ) )
      
      ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
      save( df_dday, file=filn )

    } else {
      
      load( filn )
    
    }
    if (verbose) print("done")


    ##--------------------------------------------------------
    ## re-arrange MODIS dataframe
    ##--------------------------------------------------------
    infiln  <- paste( "data/modis_nn/modis_nn_", sitename, ".Rdata", sep="" )
    outfiln <- paste( "data/df_dday_modis/df_dday_modis_", sitename, ".Rdata", sep="" )

    if (!file.exists(outfiln)||overwrite){

      if (file.exists(infiln)){

        load( infiln )  # loads 'nice_to_modis', file prepared in 'aggregate_nn_fluxnet2015.R'
        avl_modis <- TRUE
        
        droughts_modis <- get_consecutive( nice_to_modis$is_drought_byvar, leng_threshold=2, do_merge=FALSE )

        before_modis <- floor( max(droughts_modis$len) / 2 )
        after_modis  <- max(droughts_modis$len)
        
        df_dday_modis <- data.frame()
        nice_to_modis <- nice_to_modis %>% mutate( mysitename=sitename )
        for ( iinst in 1:nrow(droughts_modis) ){
          after_inst <- min( after, droughts_modis$len[iinst] )
          dday <- seq( from=-before, to=after_inst, by=1 )
          idxs <- dday + droughts_modis$idx_start[iinst]
          drophead <- which( idxs < 1 )
          if (length(drophead)>0){
            idxs <- idxs[ -drophead ]
            dday <- dday[ -drophead ]
          }
          addrows <- nice_to_modis %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst )
          df_dday_modis <- rbind( df_dday_modis, addrows )
        }

        ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
        save( df_dday_modis, file=outfiln )

        ## aggregate by 'dday'
        df_dday_aggbydday_modis <- df_dday_modis %>%  group_by( dday ) %>% 
                                                      summarise( bias_modis_med=median( bias_modis, na.rm=TRUE ), 
                                                                 bias_modis_upp=quantile( bias_modis, 0.75, na.rm=TRUE ), 
                                                                 bias_modis_low=quantile( bias_modis, 0.25, na.rm=TRUE ) ) %>%
                                                      mutate( mysitename=sitename )

        ## drop rows dday=NA
        df_dday_aggbydday_modis <- df_dday_aggbydday_modis[ which( !is.na(df_dday_aggbydday_modis$dday)), ]
        save( df_dday_aggbydday_modis, file=paste( "data/df_dday_aggbydday_modis/df_dday_aggbydday_modis_", sitename, ".Rdata", sep="" ) )

      } else {

        df_dday_modis <- NULL

      }

    } else {

      load( outfiln )

    }

    ##--------------------------------------------------------
    ## re-arrange MTE dataframe
    ##--------------------------------------------------------
    infiln  <- paste( "data/mte_nn/mte_nn_", sitename, ".Rdata", sep="" )
    outfiln <- paste( "data/df_dday_mte/df_dday_mte_", sitename, ".Rdata", sep="" )

    if (!file.exists(outfiln)||overwrite){

      if (file.exists(infiln)){

        load( infiln )
        avl_mte <- TRUE

        droughts_mte <- get_consecutive( nice_to_mte$is_drought_byvar, leng_threshold=2, do_merge=FALSE )

        if (nrow(droughts_mte)>0){

          before_mte <- floor( max(droughts_mte$len) / 2 )
          after_mte  <- max(droughts_mte$len)
    
          df_dday_mte <- data.frame()
          nice_to_mte <- nice_to_mte %>% mutate( mysitename=sitename )
          for ( iinst in 1:nrow(droughts_mte) ){
            after_inst <- min( after, droughts_mte$len[iinst] )
            dday <- seq( from=-before, to=after_inst, by=1 )
            idxs <- dday + droughts_mte$idx_start[iinst]
            drophead <- which( idxs < 1 )
            if (length(drophead)>0){
              idxs <- idxs[ -drophead ]
              dday <- dday[ -drophead ]
            }
            addrows <- nice_to_mte %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst )
            df_dday_mte <- rbind( df_dday_mte, addrows )
          }

          ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
          save( df_dday_mte, file=outfiln )

          ## aggregate by 'dday'
          df_dday_aggbydday_mte <- df_dday_mte %>%  group_by( dday ) %>% 
                                                    summarise(  bias_mte_med=median( bias_mte, na.rm=TRUE ), 
                                                                bias_mte_upp=quantile( bias_mte, 0.75, na.rm=TRUE ), 
                                                                bias_mte_low=quantile( bias_mte, 0.25, na.rm=TRUE )) %>%
                                                    mutate( mysitename=sitename )

          ## drop rows dday=NA
          df_dday_aggbydday <- df_dday_aggbydday[ which( !is.na(df_dday_aggbydday$dday)), ]
          save( df_dday_aggbydday, file=paste( "data/df_dday_aggbydday_mte/df_dday_aggbydday_mte_", sitename, ".Rdata", sep="" ) )

        } else {

          df_dday_mte <- NULL

        }

      } else {

        df_dday_mte <- NULL

      }

    } else {

      load( outfiln )

    }

  } else {

    df_dday <- NA
    df_dday_aggbydday <- NA
    df_dday_modis <- NA
    df_dday_mte <- NA

  }

  out <- list( df_dday=df_dday, df_dday_aggbydday=df_dday_aggbydday, df_dday_modis=df_dday_modis, df_dday_mte=df_dday_mte )
  return( out )

}