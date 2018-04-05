plot_fit_gpp_vs_time <- function( linearfit1, linearfit_mid, linearfit_strong, ddf=NULL, nice_agg=NULL, makepdf=FALSE ){

  require(dplyr)
  require(lubridate)

  if (!is.null(ddf)){

    ##-----------------------------------------------
    ## Time series: GPPobs and (GPP_Pmodel * fLUEest) from approaches I, II, and III
    ## aggregated to weekly
    ##-----------------------------------------------
    nice_agg <- nice_agg %>% mutate( inlowbin1 = ifelse( fvar < 0.2, 1, NA ), inlowbin2 = ifelse( fvar < 0.3, 1, NA ), inlowbin3 = ifelse( fvar < 0.4, 1, NA ) )
    wdf <- nice_agg %>% mutate( gpp_pmodel_I = gpp_pmodel * flue_est_I, 
                                gpp_pmodel_IV = gpp_pmodel * flue_est_IV , 
                                gpp_pmodel_III = gpp_pmodel * flue_est_III ) %>%
                        group_by( mysitename, week(date), year(date) ) %>% 
                        summarise( gpp_obs = mean(gpp_obs, na.rm=TRUE),
                                   gpp_pmodel  = mean(gpp_pmodel , na.rm=TRUE),
                                   gpp_pmodel_I = mean(gpp_pmodel_I, na.rm=TRUE),
                                   gpp_pmodel_IV = mean(gpp_pmodel_IV, na.rm=TRUE),
                                   gpp_pmodel_III = mean(gpp_pmodel_III, na.rm=TRUE),
                                   date = mean(date)
                                   ) %>%
                        mutate( year=year(date)) %>%
                        arrange( mysitename, year, date )

    filn <- "fig/gpp_per_site.pdf" 
    if (makepdf) print( paste( "plotting GPPobs and (GPP_Pmodel * fLUEest) vs. time for each site into file ", filn, "..." ) )
    if (makepdf) pdf( filn, width = 10, height = 6 )
    for (sitename in linearfit2$data$mysitename){
      
      df_tmp <- dplyr::filter(wdf, mysitename==sitename)
      
      if (nrow(df_tmp)>0 && any(!is.na(df_tmp$gpp_obs))){
        
        par(las=1)
        plot(  df_tmp$date, df_tmp$gpp_obs, xlab="time", ylab="GPP (gC m-2 d-1)", col=add_alpha("black", 0.5), pch=16, ylim=c( 0, max( c( df_tmp$gpp_obs, df_tmp$gpp_pmodel ), na.rm=TRUE ) ) )
        lines( df_tmp$date, df_tmp$gpp_pmodel, col="grey50" )
        lines( df_tmp$date, df_tmp$gpp_pmodel_I, col="springgreen3" )
        lines( df_tmp$date, df_tmp$gpp_pmodel_IV, col="royalblue3" )
        lines( df_tmp$date, df_tmp$gpp_pmodel_III, col="tomato" )
        title( sitename )

        legend( "topright", c("P-model", "corrected, approach I", "corrected, approach II", "corrected, approach III"), lty=1, bty="n", lwd=2, col=c("grey50", "springgreen3", "royalblue3", "tomato") )
        legend( "topleft", c("fLUE bin (0.0-0.2)", "fLUE bin (0.2-0.3)", "fLUE bin (0.3-0.4)"), pch=16, bty="n", col=c("red", "orange", "yellow") )

        with( filter(nice_agg,  mysitename==sitename), points( date, inlowbin3*0.0, pch=16, col="yellow" ) )
        with( filter(nice_agg,  mysitename==sitename), points( date, inlowbin2*0.0, pch=16, col="orange" ) )
        with( filter(nice_agg,  mysitename==sitename), points( date, inlowbin1*0.0, pch=16, col="red" ) )
      }
    }
    if (makepdf) dev.off()
  }

}