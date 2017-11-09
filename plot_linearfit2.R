plot_linearfit2 <- function( linearfit, target="fvar", df=NULL ){

  require(dplyr)

  source("stress_quad_1sided.R")

  ## merge vegetation class info into data frame
  load( "../nn_fluxnet2015/data/overview_data_fluxnet2015_L1.Rdata" ) # loads 'overview'
  linearfit$data <- linearfit$data %>% left_join( dplyr::select( overview, mysitename, classid ), by="mysitename" )

  growtype <- list( herb=c("GRA", "CRO"), sav=c("SAV", "WSA"), shrub=c("OSH", "CSH"), woody_dec=c("MF", "DBF"), woody_evg=c("ENF", "EBF"), wet=c("WET") )

  print("plotting scatter plot: fLUE0 vs. alpha (one point per site)")
  ##-----------------------------------------------
  ## Plot scatter plot: fLUE0 vs. alpha (one point per site)
  ##-----------------------------------------------
  par( las=1 )
  with( linearfit$data, plot( meanalpha, y0, pch=16, xlab="AET/PET", ylab=expression(paste("fLUE"[0])), xlim=c(0,1.1), type="n" ) )
  abline( linearfit$linmod, col="black" )

  mtext( line=-1.5, bquote( italic(R)^2 == .(format( summary( linearfit$linmod )$r.squared, digits = 2) ) ),  adj=0.1, cex=1 )
  cf <- coef(linearfit$linmod) %>% round( 2 )
  eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  mtext( line=-2.5, eq, adj=0.1 )

  ## herbaceous
  with( dplyr::filter( linearfit$data, classid %in% (growtype$herb) ), points( meanalpha, y0, pch=16, col='black', cex=1.0 ) )
  
  ## savannah
  with( dplyr::filter( linearfit$data, classid %in% (growtype$sav) ), points( meanalpha, y0, pch=18, col='black', cex=1.2 ) )
  
  ## evergreen (woody)
  with( dplyr::filter( linearfit$data, classid %in% (growtype$woody_evg) ), points( meanalpha, y0, pch=17, col='black', cex=1.0 ) )

  ## wetland
  with( dplyr::filter( linearfit$data, classid %in% (growtype$wet) ), points( meanalpha, y0, pch=25, col='black', bg='black' ) )

  ## deciduous
  with( dplyr::filter( linearfit$data, classid %in% (growtype$woody_dec) ), points( meanalpha, y0, pch=15, col='black' ) )

  ## shrublands
  with( dplyr::filter( linearfit$data, classid %in% (growtype$shrub) ), points( meanalpha, y0, pch=8, col='black' ) )

  ## label: site name
  with( linearfit$data, text( meanalpha+0.02, y0, mysitename, adj=c(0,0.5), col='black', cex=0.6 ) )

  ## legend
  legend( "bottomright", c("herbaceous", "savannah", "woody evergreen", "wetlands", "shrublands", "woody deciduous"), pch=c(16,18,17,25,8,15), bty="n", cex=0.8, inset=c(0.05,0) )


  if (!is.null(df)){
    print("plotting fLUE vs. soil moisture for each site...")
    ##-----------------------------------------------
    ## Plot scatter plot: fLUE vs. soil moisture for each site
    ##-----------------------------------------------
    for (sitename in linearfit$data$mysitename){

      df_tmp <- dplyr::filter(df, mysitename==sitename)
      data_tmp <- dplyr::filter( linearfit$data, mysitename==sitename )

      print(sitename)
      pdf( paste0( "fig/fit_", sitename, ".pdf" ) )
        par(las=1)
        plot( df_tmp$soilm_mean, df_tmp[[ target ]], xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("royalblue3", 0.2) )
        abline( h=1.0, lwd=0.5 )

        if (!is.na(dplyr::select( data_tmp, meanalpha))){ 
    
          ## Plot 1-sided curve using estimated y0 as a function of mean alpha (red)
          mycurve(  function(x) stress_quad_1sided_alpha( x, 
                                                          dplyr::select( data_tmp, meanalpha), 
                                                          dplyr::select( data_tmp, x0), 
                                                          coef(linearfit$linmod)[["(Intercept)"]], 
                                                          coef(linearfit$linmod)[["meanalpha"]] 
                                                         ),
                    from=0.0, to=1.0, col='red', add=TRUE, lwd=2 )
    
          ## Plot 1-sided curve using estimated y0 as a function of mean alpha (red)
          mycurve(  function(x) stress_quad_1sided( x, 
                                                    dplyr::select( data_tmp, x0), 
                                                    dplyr::select( data_tmp, beta)
                                                   ),
                    from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )
    
        }
      dev.off()

    }
  }

}