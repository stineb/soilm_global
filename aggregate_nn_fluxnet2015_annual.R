library(dplyr)
library(R.matlab)
library(readr)
library(lubridate)

source("get_fluxdata_fluxnet2015_annual.R")
source("remove_outliers.R")

do.siteplot <- FALSE

##------------------------------------------------
## Select all sites for which method worked (codes 1 and 2 determined by 'nn_getfail_fluxnet2015.R')
##------------------------------------------------
successcodes <- read_csv( "successcodes.csv" )

# ## Exclude sites for which no fapar data is available and hence no model results
# df_error_fapar <- read_csv( paste0( myhome, "sofun/input_fluxnet2015_sofun/error_missing_forcingdata_MODIS_FPAR_MCD15A3H_fluxnet2015.csv" ) ) 
# successcodes <- successcodes %>% left_join( df_error_fapar, by="mysitename" ) %>% rename( error_fapar = error ) %>% filter( error_fapar == 0 )

do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

##------------------------------------------------
## Load daily and 8-day data
##------------------------------------------------
filn <- paste0("data/nice_nn_agg_lue_obs_evi.Rdata")
# print( paste( "loading dataframe nice_agg from file", filn) )
load( file=filn )  # loads nice_agg
ddf <- nice_agg
rm("nice_agg")

## save aggregated NN mte
filn <- paste0("data/nice_nn_8d_agg_lue_obs_evi.Rdata")
# print( paste( "loading dataframe mte_agg from file", filn) )
load( file=filn )  # loads nice_8d_agg
ddf8 <- nice_8d_agg
rm("nice_8d_agg")

##------------------------------------------------
## Aggregate daily to annual data
##------------------------------------------------
## Probelm: Some years have very low number of modelled values. Why? Some filtering in nice_agg? check in get_modobs.R -> get_modobs_fluxnet2015.R. Change order in get_modobs.R, l.104???
number <- ddf %>%   mutate(  year = year(date) ) %>%  group_by( mysitename, year ) %>% 
                    summarise(  number_pmodel = sum(!is.na(gpp_pmodel)),
                                number_bess_v1 = sum(!is.na(gpp_bess_v1)),
                                number_bess_v2 = sum(!is.na(gpp_bess_v2))
                             )

adf <- ddf %>%  mutate(  year = year(date),
                         gpp_pmodel_s1b        = gpp_pmodel  * flue_est_IV ,
                         gpp_bess_v1_s1b       = gpp_bess_v1 * flue_est_IV ,
                         gpp_bess_v2_s1b       = gpp_bess_v2 * flue_est_IV ) %>%
                group_by( mysitename, year ) %>% 
                summarise( gpp_obs_fromdaily     = sum(gpp_obs),
                           gpp_pmodel            = sum(gpp_pmodel, na.rm=TRUE),
                           gpp_bess_v1           = sum(gpp_bess_v1),
                           gpp_bess_v2           = sum(gpp_bess_v2),
                           gpp_pmodel_s1b        = sum(gpp_pmodel_s1b, na.rm=TRUE),
                           gpp_bess_v1_s1b       = sum(gpp_bess_v1_s1b),
                           gpp_bess_v2_s1b       = sum(gpp_bess_v2_s1b)
                         ) %>%
                left_join( number, by=c("mysitename", "year") ) %>%
                mutate(  gpp_pmodel     = ifelse( number_pmodel  < 364, NA, gpp_pmodel ),
                         gpp_pmodel_s1b = ifelse( number_pmodel  < 364, NA, gpp_pmodel_s1b ),
                         gpp_bess_v1    = ifelse( number_bess_v1 < 364, NA, gpp_bess_v1 ),
                         gpp_bess_v2    = ifelse( number_bess_v2 < 364, NA, gpp_bess_v2 )
                         )

##------------------------------------------------
## Aggregate 8-daily to annual data and merge into 'adf'
##------------------------------------------------
## Probelm: Some years have very low number of modelled values. Why? Some filtering in nice_agg? check in get_modobs.R -> get_modobs_fluxnet2015.R. Change order in get_modobs.R, l.104???
number8 <- ddf8 %>%  mutate( year = year(date) ) %>%  group_by( mysitename, year ) %>% 
                     summarise( number_modis = sum(!is.na(gpp_modis)), number_vpm = sum(!is.na(gpp_vpm)) )

adf <- ddf8 %>%  mutate( year = year(date),
                          date_end = lead( date_start ),
                          gpp_modis_s1b   = gpp_modis * flue_est_IV ,
                          gpp_vpm_s1b     = gpp_vpm   * flue_est_IV ) %>%
                  mutate( ndays = ifelse( lead(mysitename)==mysitename, as.numeric( as.duration( interval( date_start, date_end ) ), "days" ), 8 ) ) %>%
                  mutate( gpp_modis     = gpp_modis     * ndays,
                          gpp_vpm       = gpp_vpm       * ndays,
                          gpp_modis_s1b = gpp_modis_s1b * ndays,
                          gpp_vpm_s1b   = gpp_vpm_s1b   * ndays
                          ) %>%
                  group_by( mysitename, year ) %>% 
                  summarise( gpp_modis     = sum( gpp_modis, na.rm=TRUE ),
                             gpp_vpm       = sum( gpp_vpm, na.rm=TRUE ),
                             gpp_modis_s1b = sum( gpp_modis_s1b, na.rm=TRUE ),
                             gpp_vpm_s1b   = sum( gpp_vpm_s1b, na.rm=TRUE )
                             ) %>%
                  left_join( number8, by = c("mysitename", "year") ) %>%
                  mutate(  gpp_modis     = ifelse( number_modis < 46, NA, gpp_modis ),
                           gpp_modis_s1b = ifelse( number_modis < 46, NA, gpp_modis_s1b ),
                           gpp_vpm       = ifelse( number_vpm < 46, NA, gpp_vpm ),
                           gpp_vpm_s1b   = ifelse( number_vpm < 46, NA, gpp_vpm_s1b )
                           ) %>%
                  right_join( adf, by = c("mysitename", "year") )

##------------------------------------------------
## Loop over sites to read annual FLUXNET data
##------------------------------------------------
## daily dataframe                  
adf_agg <- tibble()
for (sitename in unique(adf$mysitename)){

  tmp <- get_fluxdata_fluxnet2015_annual( sitename ) %>% select( year, gpp_obs ) %>% mutate( mysitename = sitename ) %>%
                                                         mutate( gpp_obs = remove_outliers( gpp_obs, coef=1.5 ) )
  addrows <- adf %>% filter( mysitename==sitename ) %>% left_join( tmp, by = c("mysitename", "year") )
  adf_agg <- bind_rows( adf_agg, addrows )

}
adf <- adf_agg
rm("adf_agg")


##------------------------------------------------
## Get statistics of annual values
##------------------------------------------------
## Calculate RMSE of pooled annual values 
stats_adf <- list(  pmodel_s0   = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_pmodel      ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_pmodel      - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    pmodel_s1b  = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_pmodel_s1b  ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_pmodel_s1b  - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    bess_v1_s0  = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_bess_v1     ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_bess_v1     - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    bess_v1_s1b = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_bess_v1_s1b ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_bess_v1_s1b - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    modis_s0    = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_modis       ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_modis       - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    vpm_s0      = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_vpm         ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_vpm         - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    modis_s1b   = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_modis_s1b   ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_modis_s1b   - adf$gpp_obs)^2, na.rm = TRUE ) ) ),
                    vpm_s1b     = list( rsq = summary( lm( adf$gpp_obs ~ adf$gpp_vpm_s1b     ) )$adj.r.squared, rmse = sqrt( mean( (adf$gpp_vpm_s1b     - adf$gpp_obs)^2, na.rm = TRUE ) ) )
                 )


##------------------------------------------------
## Get mean annual GPP per site
##------------------------------------------------
meandf <- adf %>% group_by( mysitename ) %>%
                  summarise( gpp_obs           = mean(gpp_obs, na.rm = TRUE ),
                             gpp_obs_fromdaily = mean(gpp_obs_fromdaily, na.rm = TRUE ),

                             gpp_pmodel        = mean(gpp_pmodel, na.rm = TRUE ),
                             gpp_bess_v1       = mean(gpp_bess_v1, na.rm = TRUE ),
                             gpp_bess_v2       = mean(gpp_bess_v2, na.rm = TRUE ),
                             
                             gpp_pmodel_s1b    = mean(gpp_pmodel_s1b, na.rm = TRUE ),
                             gpp_bess_v1_s1b   = mean(gpp_bess_v1_s1b, na.rm = TRUE ),
                             gpp_bess_v2_s1b   = mean(gpp_bess_v2_s1b, na.rm = TRUE ),

                             gpp_modis         = mean(gpp_modis, na.rm = TRUE ),
                             gpp_vpm           = mean(gpp_vpm  , na.rm = TRUE ),
                             
                             gpp_modis_s1b     = mean(gpp_modis_s1b, na.rm = TRUE ),
                             gpp_vpm_s1b       = mean(gpp_vpm_s1b  , na.rm = TRUE )
                  )

## Calculate RMSE of site mean values 
stats_meandf <- list(  pmodel_s0   = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_pmodel      ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_pmodel      - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       pmodel_s1b  = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_pmodel_s1b  ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_pmodel_s1b  - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       
                       bess_v1_s0  = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_bess_v1     ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_bess_v1     - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       bess_v1_s1b = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_bess_v1_s1b ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_bess_v1_s1b - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       
                       modis_s0    = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_modis       ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_modis       - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       modis_s1b   = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_modis_s1b   ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_modis_s1b   - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       
                       vpm_s0      = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_vpm         ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_vpm         - meandf$gpp_obs)^2, na.rm = TRUE ) ) ),
                       vpm_s1b     = list( rsq = summary( lm( meandf$gpp_obs ~ meandf$gpp_vpm_s1b     ) )$adj.r.squared, rmse = sqrt( mean( (meandf$gpp_vpm_s1b     - meandf$gpp_obs)^2, na.rm = TRUE ) ) )
                                            
                    )                  

##------------------------------------------------
## Get linear model for annual GPP for each site
##------------------------------------------------

linmod_list_pmodel_s0 <- list()
linmod_list_pmodel_s1b <- list()

linmod_list_bess_v1_s0 <- list()
linmod_list_bess_v1_s1b <- list()

linmod_list_modis_s0 <- list()
linmod_list_modis_s1b <- list()

linmod_list_vpm_s0 <- list()
linmod_list_vpm_s1b <- list()

slope_wgt_pmodel_s0 <- tibble()
slope_wgt_pmodel_s1b <- tibble()

slope_wgt_bess_v1_s0 <- tibble()
slope_wgt_bess_v1_s1b <- tibble()

slope_wgt_vpm_s0 <- tibble()
slope_wgt_vpm_s1b <- tibble()

slope_wgt_modis_s0 <- tibble()
slope_wgt_modis_s1b <- tibble()

for (sitename in unique(adf$mysitename)){
  
  tmp <- filter( adf, mysitename==sitename )
  
  ## pmodel
  if (sum(!is.na(tmp$gpp_pmodel))>2){

    linmod_list_pmodel_s0[[ sitename ]]  <- lm( gpp_obs ~ gpp_pmodel    , data = tmp )
    linmod_list_pmodel_s1b[[ sitename ]] <- lm( gpp_obs ~ gpp_pmodel_s1b, data = tmp )

    slope_wgt_pmodel_s0  <- slope_wgt_pmodel_s0  %>% bind_rows( tibble( slope = coef(linmod_list_pmodel_s0 [[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    slope_wgt_pmodel_s1b <- slope_wgt_pmodel_s1b %>% bind_rows( tibble( slope = coef(linmod_list_pmodel_s1b[[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    
    if (length(linmod_list_pmodel_s0[[ sitename ]]$fitted.values)>2){

      if (do.siteplot){
        par(las=1)
        with( filter( adf, mysitename==sitename ), plot(   gpp_pmodel,     gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=1, col=rgb(0,0,0,0.5) ) )
        with( filter( adf, mysitename==sitename ), points( gpp_pmodel_s1b, gpp_obs, pch=16, col=rgb(0,0,0,0.5) ) )
        lines( linmod_list_pmodel_s0[[ sitename ]]$fitted.values,  tmp$gpp_obs[ as.numeric(names(linmod_list_pmodel_s0[[ sitename ]]$fitted.values)) ],  col=rgb(1,0,0,1), lty=3 )
        lines( linmod_list_pmodel_s1b[[ sitename ]]$fitted.values, tmp$gpp_obs[ as.numeric(names(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)) ], col=rgb(1,0,0,1) )
        lines( c(0,5000), c(0,5000), lty=3 )
        mtext( sitename, line = 0.5, adj = 0, font = 2 )
        mtext( "P-model", line = 1, font = 2 )        
      }
      
    } else {

      linmod_list_pmodel_s0[[ sitename ]] <- NULL
    
    }

  }

  ## bess_v1
  if (sum(!is.na(tmp$gpp_bess_v1))>2){

    linmod_list_bess_v1_s0[[ sitename ]]  <- lm( gpp_obs ~ gpp_bess_v1, data = tmp )
    linmod_list_bess_v1_s1b[[ sitename ]] <- lm( gpp_obs ~ gpp_bess_v1_s1b, data = tmp )

    slope_wgt_bess_v1_s0  <- slope_wgt_bess_v1_s0  %>% bind_rows( tibble( slope = coef(linmod_list_bess_v1_s0 [[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    slope_wgt_bess_v1_s1b <- slope_wgt_bess_v1_s1b %>% bind_rows( tibble( slope = coef(linmod_list_bess_v1_s1b[[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    
    if (length(linmod_list_bess_v1_s0[[ sitename ]]$fitted.values)>2){

      if (do.siteplot){
        par(las=1)
        with( filter( adf, mysitename==sitename ), plot(   gpp_bess_v1, gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=1, col=rgb(0,0,0,0.5) ) )
        with( filter( adf, mysitename==sitename ), points( gpp_bess_v1_s1b, gpp_obs, pch=16, col=rgb(0,0,0,0.5) ) )
        lines( linmod_list_bess_v1_s0[[ sitename ]]$fitted.values,  tmp$gpp_obs[ as.numeric(names(linmod_list_bess_v1_s0[[ sitename ]]$fitted.values)) ],  col=rgb(1,0,0,1), lty=3 )
        lines( linmod_list_bess_v1_s1b[[ sitename ]]$fitted.values, tmp$gpp_obs[ as.numeric(names(linmod_list_bess_v1_s1b[[ sitename ]]$fitted.values)) ], col=rgb(1,0,0,1) )
        lines( c(0,5000), c(0,5000), lty=3 )
        mtext( sitename, line = 0.5, adj = 0, font = 2 )
        mtext( "BESS v1", line = 1, font = 2 )
      }
      
    } else {

      linmod_list_bess_v1_s0[[ sitename ]] <- NULL
    
    }

  }

  ## modis
  if (sum(!is.na(tmp$gpp_modis))>2){

    linmod_list_modis_s0[[ sitename ]]  <- lm( gpp_obs ~ gpp_modis    , data = tmp )
    linmod_list_modis_s1b[[ sitename ]] <- lm( gpp_obs ~ gpp_modis_s1b, data = tmp )

    slope_wgt_modis_s0  <- slope_wgt_modis_s0  %>% bind_rows( tibble( slope = coef(linmod_list_modis_s0 [[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    slope_wgt_modis_s1b <- slope_wgt_modis_s1b %>% bind_rows( tibble( slope = coef(linmod_list_modis_s1b[[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    
    if (length(linmod_list_modis_s0[[ sitename ]]$fitted.values)>2){

      if (do.siteplot){
        par(las=1)
        with( filter( adf, mysitename==sitename ), plot(   gpp_modis,     gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=1, col=rgb(0,0,0,0.5) ) )
        with( filter( adf, mysitename==sitename ), points( gpp_modis_s1b, gpp_obs, pch=16, col=rgb(0,0,0,0.5) ) )
        lines( linmod_list_modis_s0[[ sitename ]]$fitted.values,  tmp$gpp_obs[ as.numeric(names(linmod_list_modis_s0[[ sitename ]]$fitted.values)) ],  col=rgb(1,0,0,1), lty=3 )
        lines( linmod_list_modis_s1b[[ sitename ]]$fitted.values, tmp$gpp_obs[ as.numeric(names(linmod_list_modis_s1b[[ sitename ]]$fitted.values)) ], col=rgb(1,0,0,1) )
        lines( c(0,5000), c(0,5000), lty=3 )
        mtext( sitename, line = 0.5, adj = 0, font = 2 )
        mtext( "MODIS", line = 1, font = 2 )
      }
      
    } else {

      linmod_list_modis_s0[[ sitename ]] <- NULL
    
    }

  }

  ## vpm
  if (sum(!is.na(tmp$gpp_vpm))>2){

    linmod_list_vpm_s0[[ sitename ]]  <- lm( gpp_obs ~ gpp_vpm    , data = tmp )
    linmod_list_vpm_s1b[[ sitename ]] <- lm( gpp_obs ~ gpp_vpm_s1b, data = tmp )

    slope_wgt_vpm_s0  <- slope_wgt_vpm_s0  %>% bind_rows( tibble( slope = coef(linmod_list_vpm_s0 [[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    slope_wgt_vpm_s1b <- slope_wgt_vpm_s1b %>% bind_rows( tibble( slope = coef(linmod_list_vpm_s1b[[ sitename ]])[2], wgt = sd( tmp$gpp_obs ) ) )
    
    if (length(linmod_list_vpm_s0[[ sitename ]]$fitted.values)>2){

      if (do.siteplot){
        par(las=1)
        with( filter( adf, mysitename==sitename ), plot(   gpp_obs, gpp_vpm, xlim = c(0,3000), ylim = c(0,3000), pch=1, col=rgb(0,0,0,0.5) ) )
        with( filter( adf, mysitename==sitename ), points( gpp_obs, gpp_vpm_s1b, pch=16, col=rgb(0,0,0,0.5) ) )
        lines( linmod_list_vpm_s0[[ sitename ]]$fitted.values,  tmp$gpp_obs[ as.numeric(names(linmod_list_vpm_s0[[ sitename ]]$fitted.values)) ],  col=rgb(1,0,0,1), lty=3 )
        lines( linmod_list_vpm_s1b[[ sitename ]]$fitted.values, tmp$gpp_obs[ as.numeric(names(linmod_list_vpm_s1b[[ sitename ]]$fitted.values)) ], col=rgb(1,0,0,1) )
        lines( c(0,5000), c(0,5000), lty=3 )
        mtext( sitename, line = 0.5, adj = 0, font = 2 )
        mtext( "VPM", line = 1, font = 2 )
      }

    } else {

      linmod_list_vpm_s0[[ sitename ]] <- NULL
    
    }

  }

}

## Calculate mean slope weighted by SD 
stats_adf$pmodel_s0$meanslope   <- sum( slope_wgt_pmodel_s0$slope   * slope_wgt_pmodel_s0$wgt, na.rm=TRUE )   / sum( slope_wgt_pmodel_s0$wgt, na.rm=TRUE )
stats_adf$pmodel_s1b$meanslope  <- sum( slope_wgt_pmodel_s1b$slope  * slope_wgt_pmodel_s1b$wgt, na.rm=TRUE )  / sum( slope_wgt_pmodel_s1b$wgt, na.rm=TRUE )

stats_adf$bess_v1_s0$meanslope  <- sum( slope_wgt_bess_v1_s0$slope  * slope_wgt_bess_v1_s0$wgt, na.rm=TRUE )  / sum( slope_wgt_bess_v1_s0$wgt, na.rm=TRUE )
stats_adf$bess_v1_s1b$meanslope <- sum( slope_wgt_bess_v1_s1b$slope * slope_wgt_bess_v1_s1b$wgt, na.rm=TRUE ) / sum( slope_wgt_bess_v1_s1b$wgt, na.rm=TRUE )

stats_adf$modis_s0$meanslope    <- sum( slope_wgt_modis_s0$slope    * slope_wgt_modis_s0$wgt, na.rm=TRUE )    / sum( slope_wgt_modis_s0$wgt, na.rm=TRUE )
stats_adf$modis_s1b$meanslope   <- sum( slope_wgt_modis_s1b$slope   * slope_wgt_modis_s1b$wgt, na.rm=TRUE )   / sum( slope_wgt_modis_s1b$wgt, na.rm=TRUE )

stats_adf$vpm_s0$meanslope      <- sum( slope_wgt_vpm_s0$slope      * slope_wgt_vpm_s0$wgt, na.rm=TRUE )      / sum( slope_wgt_vpm_s0$wgt, na.rm=TRUE )
stats_adf$vpm_s1b$meanslope     <- sum( slope_wgt_vpm_s1b$slope     * slope_wgt_vpm_s1b$wgt, na.rm=TRUE )     / sum( slope_wgt_vpm_s1b$wgt, na.rm=TRUE )


##------------------------------------------------
## P-MODEL s0: Plot correlation of mean per site
##------------------------------------------------
plot_spatial_annual_s0 <- function( adf, meandf, linmod_list_pmodel, linmod_mean_pmodel, stats_adf, stats_meandf, slope_wgt_pmodel_s0, filn=NA ){

  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    linmod_mean_pmodel <- lm( gpp_obs ~ gpp_pmodel, data = meandf )
    stats_meandf$pmodel_s0$meanslope <- coef(linmod_mean_pmodel)[2]
    par(las=1, mar=c(4,4.5,4,1))

    # with( adf, plot(   gpp_obs, gpp_pmodel, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
    # abline( linmod_mean_pmodel, col="red")
    # lines( c(0,5000), c(0,5000), lty=3 )

    with( meandf, plot( gpp_pmodel, gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", ylab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), xlab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
    abline( linmod_mean_pmodel, col="red")
    for (sitename in ls(linmod_list_pmodel) ){
      tmp <- filter( adf, mysitename==sitename )
      if (length(linmod_list_pmodel[[ sitename ]]$fitted.values)>3){
        par(las=1)
        lines( tmp$gpp_pmodel[ as.numeric(names(linmod_list_pmodel[[ sitename ]]$fitted.values)) ], linmod_list_pmodel[[ sitename ]]$fitted.values, col="black")
      }
    }
    lines( c(0,5000), c(0,5000), lty=3 )
    mtext( "P-model, s0", line = 1, font = 2 )

    bquote( italic(R)^2 == .(format( stats_adf$pmodel_s0$rsq, digits = 2) ) )

    mtext( bquote( italic(R)^2 == .(format( stats_adf$pmodel_s0$rsq, digits = 2) ) ), adj = 1, cex = 0.8, line=2 )
    mtext( paste0( "RMSE = ",  format( stats_adf$pmodel_s0$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
    # mtext( paste0( "slope = ", format( stats_adf$pmodel_s0$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

    mtext( bquote( italic(R)^2 == .(format( stats_meandf$pmodel_s0$rsq, digits = 2) ) ), adj = 0, cex = 0.8, line=2, col="red" )
    mtext( paste0( "RMSE = ",  format( stats_meandf$pmodel_s0$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
    mtext( paste0( "slope = ", format( stats_meandf$pmodel_s0$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )

    ##-----------------------------------------------------
    ## Inset: histogram of slopes
    ##-----------------------------------------------------
    ## Inset 1
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v_orig <- v
    v <- c( v[1]+0.03, v[1]+0.2*v[2], v[3]+0.50*v[4], v[3]+0.72*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.5,0) )
    hist( slope_wgt_pmodel_s0$slope, xlim=c(-2,3), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 150, xlab="slope"  )
    abline( v=1.0, col="red" )
    axis( 1, cex.axis=0.7, xlab="slope" )


  if (!is.na(filn)) dev.off()


}

plot_spatial_annual_s0( adf, meandf, linmod_list_pmodel_s0,  linmod_mean_pmodel_s0,  stats_adf, stats_meandf, slope_wgt_pmodel_s0, filn="fig/modobs_ann_spatial_temporal_pmodel_s0.pdf"  )


##------------------------------------------------
## P-MODEL s1a: Plot correlation of mean per site
##------------------------------------------------
plot_spatial_annual_s1 <- function( adf, meandf, linmod_list_pmodel, linmod_mean_pmodel, stats_adf, stats_meandf, slope_wgt_pmodel_s1b, filn=NA ){

  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    linmod_mean_pmodel_s1b <- lm( gpp_obs ~ gpp_pmodel_s1b, data = meandf )
    stats_meandf$pmodel_s1b$meanslope <- coef(linmod_mean_pmodel_s1b)[2]
    par(las=1, mar=c(4,4.5,4,1))

    # with( adf, plot(   gpp_obs, gpp_pmodel_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
    # abline( linmod_mean_pmodel_s1b, col="red")
    # lines( c(0,5000), c(0,5000), lty=3 )

    with( meandf, plot( gpp_pmodel_s1b, gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", ylab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), xlab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
    abline( linmod_mean_pmodel_s1b, col="red")
    for (sitename in ls(linmod_list_pmodel_s1b) ){
      tmp <- filter( adf, mysitename==sitename )
      if (length(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)>3){
        par(las=1)
        lines( tmp$gpp_pmodel_s1b[ as.numeric(names(linmod_list_pmodel[[ sitename ]]$fitted.values)) ], linmod_list_pmodel[[ sitename ]]$fitted.values, col="black")
      }
    }
    lines( c(0,5000), c(0,5000), lty=3 )
    mtext( "P-model, s1b", line = 1, font = 2 )
    mtext( bquote( italic(R)^2 == .(format( stats_adf$pmodel_s1b$rsq, digits = 2) ) ), adj = 1, cex = 0.8, line=2 )
    mtext( paste0( "RMSE = ", format( stats_adf$pmodel_s1b$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
    # mtext( paste0( "slope = ", format( stats_adf$pmodel_s1b$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

    mtext( bquote( italic(R)^2 == .(format( stats_meandf$pmodel_s1b$rsq, digits = 2) ) ), adj = 0, cex = 0.8, line=2, col="red" )
    mtext( paste0( "RMSE = ", format( stats_meandf$pmodel_s1b$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
    mtext( paste0( "slope = ", format( stats_meandf$pmodel_s1b$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
    
    ##-----------------------------------------------------
    ## Inset: histogram of slopes
    ##-----------------------------------------------------
    ## Inset 1
    u <- par("usr")
    v <- c(
      grconvertX(u[1:2], "user", "ndc"),
      grconvertY(u[3:4], "user", "ndc")
    )
    v_orig <- v
    v <- c( v[1]+0.03, v[1]+0.2*v[2], v[3]+0.50*v[4], v[3]+0.72*v[4] )
    par( fig=v, new=TRUE, mar=c(0,0,0,0), mgp=c(3,0.5,0) )
    hist( slope_wgt_pmodel_s1b$slope, xlim=c(-2,3), cex.axis=0.7, axes=FALSE, col="grey70", main="", breaks = 50, xlab="slope" )
    abline( v=1.0, col="red" )
    axis( 1, cex.axis=0.7, xlab="slope" )
    
  if (!is.na(filn)) dev.off()

}

plot_spatial_annual_s1( adf, meandf, linmod_list_pmodel_s1b,  linmod_mean_pmodel_s1b,  stats_adf, stats_meandf, slope_wgt_pmodel_s1b, filn="fig/modobs_ann_spatial_temporal_pmodel_s1b.pdf"  )


# ## for each site:
# linmod_mean_pmodel_s1b <- lm( gpp_obs ~ gpp_pmodel_s1b, data = meandf )
# stats_meandf$pmodel_s1b$meanslope <- coef(linmod_mean_pmodel_s1b)[2]
# par(las=1, mar=c(4,4,4,1))

# for (sitename in ls(linmod_list_pmodel_s1b) ){
#   with( meandf, plot( gpp_pmodel_s1b, gpp_obs, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", ylab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), xlab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_pmodel_s1b, col="red")
#   tmp <- filter( adf, mysitename==sitename )
#   if (length(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)>3){
#     par(las=1)
#     lines(  tmp$gpp_pmodel[ as.numeric(names(linmod_list_pmodel_s0[[ sitename ]]$fitted.values)) ], linmod_list_pmodel_s0[[ sitename ]]$fitted.values, col="grey70")
#     points( tmp$gpp_pmodel[ as.numeric(names(linmod_list_pmodel_s0[[ sitename ]]$fitted.values)) ], tmp$gpp_obs[ as.numeric(names(linmod_list_pmodel_s0[[ sitename ]]$fitted.values)) ], col="grey70", pch=16)

#     lines(  tmp$gpp_pmodel_s1b[ as.numeric(names(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)) ], linmod_list_pmodel_s1b[[ sitename ]]$fitted.values, col="black")
#     points( tmp$gpp_pmodel_s1b[ as.numeric(names(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)) ], tmp$gpp_obs[ as.numeric(names(linmod_list_pmodel_s1b[[ sitename ]]$fitted.values)) ], col="black", pch=16)
#   }
#   title( sitename )
# }



# ##------------------------------------------------
# ## BESS s0: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_bess_s0.pdf", width = 8, height = 6)
#   linmod_mean_bess_v1_s0 <- lm( gpp_bess_v1 ~ gpp_obs, data = meandf )
#   stats_meandf$bess_v1_s0$meanslope <- coef(linmod_mean_bess_v1_s0)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_bess_v1, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_bess_v1_s0, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_bess_v1, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_bess_v1_s0, col="red")
#   for (sitename in ls(linmod_list_bess_v1_s0) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_bess_v1_s0[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_bess_v1_s0[[ sitename ]]$fitted.values)) ], linmod_list_bess_v1_s0[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "BESS, s0", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$bess_v1_s0$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$bess_v1_s0$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$bess_v1_s0$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$bess_v1_s0$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()

# ##------------------------------------------------
# ## BESS s1a: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_bess_s1b.pdf", width = 8, height = 6)
#   linmod_mean_bess_v1_s1b <- lm( gpp_bess_v1_s1b ~ gpp_obs, data = meandf )
#   stats_meandf$bess_v1_s1b$meanslope <- coef(linmod_mean_bess_v1_s1b)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_bess_v1_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_bess_v1_s1b, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_bess_v1_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_bess_v1_s1b, col="red")
#   for (sitename in ls(linmod_list_bess_v1_s1b) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_bess_v1_s1b[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_bess_v1_s1b[[ sitename ]]$fitted.values)) ], linmod_list_bess_v1_s1b[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "BESS, s1b", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$bess_v1_s1b$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$bess_v1_s1b$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$bess_v1_s1b$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$bess_v1_s1b$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()

# ##------------------------------------------------
# ## MODIS s0: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_modis_s0.pdf", width = 8, height = 6)
#   linmod_mean_modis_s0 <- lm( gpp_modis ~ gpp_obs, data = meandf )
#   stats_meandf$modis_s0$meanslope <- coef(linmod_mean_modis_s0)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_modis, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_modis_s0, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_modis, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_modis_s0, col="red")
#   for (sitename in ls(linmod_list_modis_s0) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_modis_s0[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_modis_s0[[ sitename ]]$fitted.values)) ], linmod_list_modis_s0[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "MODIS, s0", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$modis_s0$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$modis_s0$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$modis_s0$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$modis_s0$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()

# ##------------------------------------------------
# ## MODIS s1a: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_modis_s1b.pdf", width = 8, height = 6)
#   linmod_mean_modis_s1b <- lm( gpp_modis_s1b ~ gpp_obs, data = meandf )
#   stats_meandf$modis_s1b$meanslope <- coef(linmod_mean_modis_s1b)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_modis_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_modis_s1b, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_modis_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_modis_s1b, col="red")
#   for (sitename in ls(linmod_list_modis_s1b) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_modis_s1b[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_modis_s1b[[ sitename ]]$fitted.values)) ], linmod_list_modis_s1b[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "MODIS, s1b", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$modis_s1b$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$modis_s1b$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$modis_s1b$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$modis_s1b$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()

# ##------------------------------------------------
# ## VPM s0: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_vpm_s0.pdf", width = 8, height = 6)
#   linmod_mean_vpm_s0 <- lm( gpp_vpm ~ gpp_obs, data = meandf )
#   stats_meandf$vpm_s0$meanslope <- coef(linmod_mean_vpm_s0)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_vpm, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_vpm_s0, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_vpm, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_vpm_s0, col="red")
#   for (sitename in ls(linmod_list_vpm_s0) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_vpm_s0[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_vpm_s0[[ sitename ]]$fitted.values)) ], linmod_list_vpm_s0[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "VPM, s0", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$vpm_s0$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$vpm_s0$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$vpm_s0$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$vpm_s0$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()

# ##------------------------------------------------
# ## VPM s1b: Plot correlation of mean per site
# ##------------------------------------------------
# pdf("fig/modobs_ann_spatial_temporal_vpm_s1b.pdf", width = 8, height = 6)
#   linmod_mean_vpm_s1b <- lm( gpp_vpm_s1b ~ gpp_obs, data = meandf )
#   stats_meandf$vpm_s1b$meanslope <- coef(linmod_mean_vpm_s1b)[2]
#   par(las=1)

#   # with( adf, plot(   gpp_obs, gpp_vpm_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5) ) )
#   # abline( linmod_mean_vpm_s1b, col="red")
#   # lines( c(0,5000), c(0,5000), lty=3 )

#   with( meandf, plot( gpp_obs, gpp_vpm_s1b, xlim = c(0,3000), ylim = c(0,3000), pch=16, col=rgb(0,0,0,0.5), type = "n", xlab = expression( paste("observed GPP (gC m"^-2, "yr"^-1, ")" ) ), ylab = expression( paste("simulated GPP (gC m"^-2, "yr"^-1, ")" ) ) ) )
#   abline( linmod_mean_vpm_s1b, col="red")
#   for (sitename in ls(linmod_list_vpm_s1b) ){
#     tmp <- filter( adf, mysitename==sitename )
#     if (length(linmod_list_vpm_s1b[[ sitename ]]$fitted.values)>3){
#       par(las=1)
#       lines( tmp$gpp_obs[ as.numeric(names(linmod_list_vpm_s1b[[ sitename ]]$fitted.values)) ], linmod_list_vpm_s1b[[ sitename ]]$fitted.values, col="black")
#     }
#   }
#   lines( c(0,5000), c(0,5000), lty=3 )
#   mtext( "VPM, s1b", line = 1, font = 2 )
#   mtext( paste0( "RMSE = ", format( stats_adf$vpm_s1b$rmse, digits = 3 ) ), adj = 1, cex = 0.8, line=1 )
#   mtext( paste0( "slope = ", format( stats_adf$vpm_s1b$meanslope, digits = 3 ) ), adj = 1, cex = 0.8 )

#   mtext( paste0( "RMSE = ", format( stats_meandf$vpm_s1b$rmse, digits = 3 ) ), adj = 0, cex = 0.8, line=1, col="red" )
#   mtext( paste0( "slope = ", format( stats_meandf$vpm_s1b$meanslope, digits = 3 ) ), adj = 0, cex = 0.8, col="red" )
# dev.off()




## define some IAV goodness metric: mean slope, weighted by standard deviation in gpp_obs; complementary information to RMSE!
