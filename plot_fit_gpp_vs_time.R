plot_fit_gpp_vs_time <- function( linearfit1, linearfit2, linearfit3, ddf=NULL, nice_agg=NULL, makepdf=FALSE ){

  require(dplyr)
  require(lubridate)

  if (!is.null(ddf)){

    ##-----------------------------------------------
    ## Time series: GPPobs and (GPP_Pmodel * fLUEest) from approaches I, II, and III
    ## aggregated to weekly
    ##-----------------------------------------------
    wdf <- nice_agg %>% mutate( gpp_pmodel1 = gpp_pmodel * flue_est_1, 
                                gpp_pmodel2 = gpp_pmodel * flue_est_2, 
                                gpp_pmodel3 = gpp_pmodel * flue_est_3 ) %>%
                        group_by( mysitename, week(date), year(date) ) %>% 
                        summarise( gpp_obs = mean(gpp_obs, na.rm=TRUE),
                                   gpp_pmodel  = mean(gpp_pmodel , na.rm=TRUE),
                                   gpp_pmodel1 = mean(gpp_pmodel1, na.rm=TRUE),
                                   gpp_pmodel2 = mean(gpp_pmodel2, na.rm=TRUE),
                                   gpp_pmodel3 = mean(gpp_pmodel3, na.rm=TRUE),
                                   date = mean(date)
                                   ) %>%
                        mutate( year=year(date)) %>%
                        arrange( mysitename, year, date )

    filn <- "fig/gpp_per_site.pdf" 
    if (makepdf) print( paste( "plotting GPPobs and (GPP_Pmodel * fLUEest) vs. time for each site into file ", filn, "..." ) )
    if (makepdf) pdf( filn, width = 10, height = 6 )
    for (sitename in linearfit2$data$mysitename){

      df_tmp <- dplyr::filter(wdf, mysitename==sitename)

      if (nrow(df_tmp)>0){

        par(las=1)
        plot(  df_tmp$date, df_tmp$gpp_obs, xlab="time", ylab="GPP (gC m-2 d-1)", col=add_alpha("black", 0.5), pch=16 )
        lines( df_tmp$date, df_tmp$gpp_pmodel, col="grey50" )
        lines( df_tmp$date, df_tmp$gpp_pmodel1, col="springgreen3" )
        lines( df_tmp$date, df_tmp$gpp_pmodel2, col="royalblue3" )
        lines( df_tmp$date, df_tmp$gpp_pmodel3, col="tomato" )
        title( sitename )
        legend( "topright", c("P-model", "corrected, approach I", "corrected, approach II", "corrected, approach III"), lty=1, bty="n", lwd=2, col=c("grey50", "springgreen3", "royalblue3", "tomato") )
      }
    }
    if (makepdf) dev.off()
  }

}