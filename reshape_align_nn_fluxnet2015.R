reshape_align_nn_fluxnet2015 <- function( sitename, nam_target="lue_obs_evi", bysm=FALSE, use_fapar=FALSE, use_weights=FALSE, overwrite=TRUE, verbose=FALSE ){

  # ## debug-------------------
  # sitename = "IT-Ro1"
  # nam_target="lue_obs_evi"
  # bysm=FALSE
  # use_fapar=FALSE
  # use_weights=FALSE
  # overwrite=TRUE
  # verbose=TRUE
  # #--------------------------

  syshome <- Sys.getenv( "HOME" )
  source( paste( syshome, "/.Rprofile", sep="" ) )

  require( dplyr )
  require( tidyr )

  source( paste( myhome, "sofun/utils_sofun/analysis_sofun/get_consecutive.R", sep="" ) ) 

  before <- 30
  after  <- 300

  ## Bins for different variables
  fvarbins  <- seq( from=-20, to=40, by=20 )
  faparbins <- seq( from=-20, to=40, by=20 )
  iwuebins  <- seq( from=-30, to=60, by=30 )

  bincentres_fvar  <- fvarbins[1:(length(fvarbins)-1)]   + (fvarbins[2]-fvarbins[1])/2
  bincentres_fapar <- faparbins[1:(length(faparbins)-1)] + (faparbins[2]-faparbins[1])/2
  bincentres_iwue  <- iwuebins[1:(length(iwuebins)-1)]   + (iwuebins[2]-iwuebins[1])/2

  usecols <- c( "mysitename",
                "date",
                "is_drought_byvar", 
                "vpd", 
                "soilm_mean", 
                "soilm_splash220", 
                "fpar", 
                "fvar", 
                "gpp_obs",
                "flue_est_1",
                "flue_est_2",
                "flue_est_3",
                "alpha", 
                "gpp_pmodel",
                "bias_pmodel", 
                "ratio_obs_mod_pmodel",
                "gpp_bess_v1", "gpp_bess_v2",
                "bias_bess_v1", "bias_bess_v2",
                "ratio_obs_mod_bess_v1", "ratio_obs_mod_bess_v2"
              )

  usecols_8d <- c(  "mysitename",
                    "date",
                    "is_drought_byvar", 
                    "vpd", 
                    "soilm_mean", 
                    "soilm_splash220", 
                    "fpar", 
                    "fvar", 
                    "gpp_obs",
                    "flue_est_1",
                    "flue_est_2",
                    "flue_est_3",
                    "alpha", 
                    "gpp_modis", 
                    "bias_modis", 
                    "ratio_obs_mod_modis",
                    "gpp_vpm",
                    "bias_vpm",
                    "ratio_obs_mod_vpm", 
                    "gpp_mte", 
                    "bias_mte", 
                    "ratio_obs_mod_mte",
                    "gpp_pmodel", 
                    "bias_pmodel", 
                    "ratio_obs_mod_pmodel",
                    "gpp_bess_v1", "gpp_bess_v2",
                    "bias_bess_v1", "bias_bess_v2",
                    "ratio_obs_mod_bess_v1", "ratio_obs_mod_bess_v2"
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
  if (!dir.exists("./data/droughts")) system( "mkdir -p ./data/droughts" )
  if (!file.exists(droughtfil)){
    dir <- paste( myhome, "data/nn_fluxnet/fvar/", sep="" )
    infil <- paste( dir, "nn_fluxnet2015_", sitename, "_lue_obs_evi.Rdata", sep="" ) 
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
    if (!dir.exists("./data/df_dday")) system( "mkdir -p ./data/df_dday" )
    
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
        addrows <- df %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst ) %>% select( mysitename, dday, inst, one_of(usecols) )
        df_dday <- df_dday %>% bind_rows( addrows )              
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
                         summarise( fpar  = median( fpar , na.rm=TRUE ) ) %>%
                         complete( infaparbin, fill = list( fpar  = NA ) ) %>% 
                         dplyr::select( fpar )
      tmp <- unlist( tmp )[1:(length(faparbins)-1)]
      df_dday$fpar_norm = df_dday$fpar / tmp[1]


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

                                                  ## FPAR
                                                  fpar_med=median( fpar, na.rm=TRUE ), fpar_upp=quantile( fpar, 0.75, na.rm=TRUE ), fpar_low=quantile( fpar, 0.25, na.rm=TRUE ),

                                                  ## P-model bias
                                                  bias_pmodel_med=median( bias_pmodel, na.rm=TRUE ), bias_pmodel_upp=quantile( bias_pmodel, 0.75, na.rm=TRUE ), bias_pmodel_low=quantile( bias_pmodel, 0.25, na.rm=TRUE ),

                                                  ## BESS v1 bias
                                                  bias_bess_v1_med=median( bias_bess_v1, na.rm=TRUE ), bias_bess_v1_upp=quantile( bias_bess_v1, 0.75, na.rm=TRUE ), bias_bess_v1_low=quantile( bias_bess_v1, 0.25, na.rm=TRUE ),

                                                  ## BESS v2 bias
                                                  bias_bess_v2_med=median( bias_bess_v2, na.rm=TRUE ), bias_bess_v2_upp=quantile( bias_bess_v2, 0.75, na.rm=TRUE ), bias_bess_v2_low=quantile( bias_bess_v2, 0.25, na.rm=TRUE )

                                                  ) %>%
                                        mutate( mysitename=sitename )

      ## drop rows dday=NA
      df_dday_aggbydday <- df_dday_aggbydday[ which( !is.na(df_dday_aggbydday$dday)), ]
      if (!dir.exists("data/df_dday_aggbydday")) system( "mkdir -p data/df_dday_aggbydday")
      save( df_dday_aggbydday, file=paste( "data/df_dday_aggbydday/df_dday_aggbydday_", sitename, ".Rdata", sep="" ) )
      
      ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
      save( df_dday, file=filn )

    } else {
      
      load( filn )
    
    }
    if (verbose) print("done")


    ##--------------------------------------------------------
    ## re-arrange 8d dataframe
    ##--------------------------------------------------------
    infiln  <- paste0( "data/nice_8d_nn/nice_8d_nn", sitename, ".Rdata"  )
    outfiln <- paste0( "data/df_dday_8d/df_dday_8d_", sitename, ".Rdata", sep="" )
    if (!dir.exists("./data/df_dday_8d")) system( "mkdir -p ./data/df_dday_8d" )
    
    if (!file.exists(outfiln)||overwrite){

      if (file.exists(infiln)){

        load( infiln )  # loads 'nice_8d', file prepared in 'aggregate_nn_fluxnet2015.R'
        avl_modis <- TRUE
        
        droughts_8d <- get_consecutive( nice_8d$is_drought_byvar, leng_threshold=2, do_merge=FALSE )
        
        df_dday_8d <- data.frame()
        nice_8d <- nice_8d %>% mutate( mysitename=sitename )
        for ( iinst in 1:nrow(droughts_8d) ){
          after_inst <- min( after, droughts_8d$len[iinst] )
          dday <- seq( from=-before, to=after_inst, by=1 )
          idxs <- dday + droughts_8d$idx_start[iinst]
          drophead <- which( idxs < 1 )
          if (length(drophead)>0){
            idxs <- idxs[ -drophead ]
            dday <- dday[ -drophead ]
          }
          addrows <- nice_8d %>% slice( idxs ) %>% mutate( dday=dday, inst=iinst ) %>% select( mysitename, dday, inst, one_of(usecols_8d ) )
          df_dday_8d <- df_dday_8d %>% bind_rows( addrows )
        }

        ## Append to Rdata file that already has the aligned array. Function 'resave()' is in my .Rprofile
        save( df_dday_8d, file=outfiln )

        ## aggregate by 'dday'
        df_dday_8d_agg_med <- df_dday_8d %>%  group_by( dday ) %>% 
                                              summarise_at( vars( starts_with("bias_") ), funs(  median(.) ), na.rm = TRUE ) %>%
                                              setNames( c("dday", paste0( names(.)[-1], "_med" ) ) )

        df_dday_8d_agg_q25 <- df_dday_8d %>%  group_by( dday ) %>% 
                                              summarise_at( vars( starts_with("bias_") ), funs(  quantile(., probs=0.25) ), na.rm = TRUE ) %>%
                                              setNames( c("dday", paste0( names(.)[-1], "_q25" ) ) )

        df_dday_8d_agg_q75 <- df_dday_8d %>%  group_by( dday ) %>% 
                                              summarise_at( vars( starts_with("bias_") ), funs(  quantile(., probs=0.75) ), na.rm = TRUE ) %>%
                                              setNames( c("dday", paste0( names(.)[-1], "_q75" ) ) )

        df_dday_aggbydday_8d <- df_dday_8d_agg_med %>%  left_join( df_dday_8d_agg_q25, by = "dday" ) %>% left_join( df_dday_8d_agg_q75, by = "dday" ) %>%
                                                        mutate( mysitename=sitename ) %>%
                                                        filter( !is.na( dday ) )

        if (!dir.exists("data/df_dday_aggbydday_8d")) system( "mkdir -p data/df_dday_aggbydday_8d")
        save( df_dday_aggbydday_8d, file=paste( "data/df_dday_aggbydday_8d/df_dday_aggbydday_8d_", sitename, ".Rdata", sep="" ) )

      } else {

        df_dday_8d <- NULL

      }

    } else {

      load( outfiln )

    }

  } else {

    df_dday           <- NULL
    df_dday_aggbydday <- NULL
    df_dday_8d        <- NULL

  }

  out <- list( df_dday=df_dday, df_dday_aggbydday=df_dday_aggbydday, df_dday_8d=df_dday_8d )
  return( out )

}

  