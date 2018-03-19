plot_linearfit3 <- function( linearfit1, linearfit2, linearfit3, ddf=NULL, nice_agg=NULL ){

  require(dplyr)
  require(lubridate)

  # ## merge vegetation class info into data frame
  # load( "../nn_fluxnet2015/data/overview_data_fluxnet2015_L1.Rdata" ) # loads 'overview'
  # linearfit$data <- linearfit$data %>% left_join( dplyr::select( overview, mysitename, classid ), by="mysitename" )

  # growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

  # ##-----------------------------------------------
  # ## Plot scatter plot: fLUE0 vs. alpha (one point per site)
  # ##-----------------------------------------------
  # par( las=1 )
  # with( linearfit$data, plot( meanalpha, y0, pch=16, xlab="AET/PET", ylab=expression(paste("fLUE"[0])), xlim=c(0,1.1), type="n" ) )
  # abline( linearfit$linmod, col="black", lty=2 )
  # abline( linearfit2$linmod, col="black" )

  # mtext( line=-1.5, bquote( italic(R)^2 == .(format( summary( linearfit2$linmod )$r.squared, digits = 2) ) ),  adj=0.1, cex=1 )
  # cf <- coef(linearfit2$linmod) %>% round( 2 )
  # eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  # mtext( line=-2.5, eq, adj=0.1 )

  # ## herbaceous
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$herb) ), points( meanalpha, y0, pch=16, col='black', cex=1.0 ) )
  
  # ## savannah
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$sav) ), points( meanalpha, y0, pch=18, col='black', cex=1.2 ) )
  
  # ## evergreen (woody)
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$woody_evg) ), points( meanalpha, y0, pch=17, col='black', cex=1.0 ) )

  # ## wetland
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$wet) ), points( meanalpha, y0, pch=25, col='black', bg='black' ) )

  # ## deciduous
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$woody_dec) ), points( meanalpha, y0, pch=15, col='black' ) )

  # ## shrublands
  # with( dplyr::filter( linearfit$data, classid %in% (growtype$shrub) ), points( meanalpha, y0, pch=8, col='black' ) )

  # ## label: site name
  # with( linearfit$data, text( meanalpha+0.02, y0, mysitename, adj=c(0,0.5), col='black', cex=0.6 ) )

  # ## legend
  # legend( "bottomright", c("herbaceous", "savannah", "woody evergreen", "wetlands", "shrublands", "woody deciduous"), pch=c(16,18,17,25,8,15), bty="n", cex=0.8, inset=c(0.05,0) )


  if (!is.null(ddf)){

    ##-----------------------------------------------
    ## Plot scatter plot: fLUE vs. soil moisture for each site
    ##-----------------------------------------------
    filn <- "fig/fit_to_bias_plot_per_site.pdf" 
    print( paste( "plotting fLUE vs. soil moisture for each site into file ", filn, "..." ) )
    pdf( filn, width = 5, height = 4 )
    for (sitename in linearfit2$data$mysitename){

      df_tmp <- dplyr::filter(ddf, mysitename==sitename)
      data_tmp <- dplyr::filter( linearfit2$data, mysitename==sitename )

      par(las=1)
      plot( df_tmp$soilm_mean, df_tmp[[ "fvar" ]], xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("black", 0.2) )
      abline( h=1.0, lwd=0.5 )
      title( sitename )

      if (!is.na(dplyr::select( data_tmp, meanalpha))){ 
  
        ## Curve from approach 1
        mycurve(  function(x) calc_flue_est_alpha(  x, 
                                                    alpha=dplyr::select( data_tmp, meanalpha), 
                                                    apar=coef(linearfit1$linmod)[1], 
                                                    bpar=coef(linearfit1$linmod)[2], 
                                                    cpar=0.125, 
                                                    dpar=0.75 
                                                    ),
                  from=0.0, to=1.0, col='springgreen3', add=TRUE, lwd=2 )

        ## Curve from approach 2
        mycurve(  function(x) stress_quad_1sided_alpha( x, 
                                                        dplyr::select( data_tmp, meanalpha), 
                                                        x0=0.9, 
                                                        coef(linearfit2$linmod)[["(Intercept)"]], 
                                                        coef(linearfit2$linmod)[["meanalpha"]] 
                                                       ),
                  from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )

        ## Curve from approach 3
        mycurve(  function(x) stress_quad_1sided_alpha( x, 
                                                        dplyr::select( data_tmp, meanalpha), 
                                                        x0=0.9, 
                                                        coef(linearfit3$linmod)[["(Intercept)"]], 
                                                        coef(linearfit3$linmod)[["meanalpha"]] 
                                                        ),
                  from=0.0, to=1.0, col='tomato', add=TRUE, lwd=2 )
        # mycurve(  function(x) stress_quad_1sided_alpha( x, 
        #                                                 dplyr::select( data_tmp, meanalpha), 
        #                                                 x0=0.9, 
        #                                                 coef(linearfit3)[["apar"]], 
        #                                                 coef(linearfit3)[["bpar"]] 
        #                                                ),
        #           from=0.0, to=1.0, col='tomato', add=TRUE, lwd=2 )
        legend( "bottomright", c("approach I", "approach II", "approach III"), lty=1, bty="n", lwd=2, col=c("springgreen3", "royalblue3", "tomato") )
  
      }
    }
    dev.off()

    ##-----------------------------------------------
    ## Time series: fLUE and fLUEest from approaches I, II, and III
    ##-----------------------------------------------
    filn <- "fig/flue_est_per_site.pdf" 
    print( paste( "plotting fLUE and fLUEest vs. time for each site into file ", filn, "..." ) )
    pdf( filn, width = 10, height = 6 )
    for (sitename in linearfit2$data$mysitename){

      df_tmp <- dplyr::filter(nice_agg, mysitename==sitename)

      if (nrow(df_tmp)>0){
        
        par(las=1)
        plot(  df_tmp$date, df_tmp[[ "fvar" ]], type="l", xlab="time", ylab="fLUE", col="black", ylim = c(0.2, 1.2) )
        lines( df_tmp$date, df_tmp$flue_est_1, col="springgreen3" )
        lines( df_tmp$date, df_tmp$flue_est_2, col="royalblue3" )
        lines( df_tmp$date, df_tmp$flue_est_3, col="tomato" )
        title( sitename )
        legend( "bottomright", c("approach I", "approach II", "approach III"), lty=1, bty="n", lwd=2, col=c("springgreen3", "royalblue3", "tomato") )
      }
    }
    dev.off()

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
    print( paste( "plotting GPPobs and (GPP_Pmodel * fLUEest) vs. time for each site into file ", filn, "..." ) )
    pdf( filn, width = 10, height = 6 )
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
    dev.off()
  }

}