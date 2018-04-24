plot_fit_vs_soilmoist <- function( linearfit_low, linearfit_mid, linearfit_strong, ddf=NULL, nice_agg=NULL, makepdf=FALSE ){

  require(dplyr, quietly = TRUE, warn.conflicts = FALSE )
  require(lubridate, quietly = TRUE, warn.conflicts = FALSE )

  source("stress_exp.R")
  source("calc_flue_est_alpha.R")
  source("stress_quad_1sided.R")

  if (!is.null(ddf)){

    ##-----------------------------------------------
    ## Plot scatter plot: fLUE vs. soil moisture for each site
    ##-----------------------------------------------
    filn <- "fig/fit_to_bias_plot_per_site.pdf" 
    if (makepdf) print( paste( "plotting fLUE vs. soil moisture for each site into file ", filn, "..." ) )
    if (makepdf) pdf( filn, width = 5, height = 4 )
    for (sitename in linearfit_mid$data$mysitename){

      df_tmp <- dplyr::filter(ddf, mysitename==sitename)
      data_tmp <- dplyr::filter( linearfit_mid$data, mysitename==sitename )

      classid <- unique(df_tmp$classid)

      par(las=1)
      plot( df_tmp$soilm_mean, df_tmp[[ "fvar" ]], xlim=c(0,1), ylim=c(0,1.2), pch=16, xlab="soil water content (fraction)", ylab="fLUE", col=add_alpha("black", 0.2) )
      abline( h=1.0, lwd=0.5 )
      mtext( sitename, line = 0.25, font = 2, adj = 0 )

      if (!is.na(dplyr::select( data_tmp, meanalpha))){ 
  
        ## Curve from approach I
        mycurve(  function(x) calc_flue_est_alpha(  x, 
                                                    alpha=dplyr::select( data_tmp, meanalpha), 
                                                    apar=coef(linearfit_low$linmod)[1], 
                                                    bpar=coef(linearfit_low$linmod)[2], 
                                                    cpar=0.125, 
                                                    dpar=0.75 
                                                    ),
                  from=0.0, to=1.0, col='springgreen3', add=TRUE, lwd=2 )

        ## Curve from approach IV (mid)
        mycurve(  function(x) stress_quad_1sided_alpha_grasstree( x, 
                                                        dplyr::select( data_tmp, meanalpha), 
                                                        x0=0.9, 
                                                        c(coef(linearfit_mid$linmod_tree)[["(Intercept)"]],coef(linearfit_mid$linmod_grass)[["(Intercept)"]]), 
                                                        c(coef(linearfit_mid$linmod_tree)[["meanalpha"]],coef(linearfit_mid$linmod_grass)[["meanalpha"]]),
                                                        classid
                                                       ),
                  from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )

        ## Curve from approach III (strong)
        mycurve(  function(x) stress_quad_1sided_alpha( x, 
                                                        dplyr::select( data_tmp, meanalpha), 
                                                        x0=0.9, 
                                                        coef(linearfit_strong$linmod)[["(Intercept)"]], 
                                                        coef(linearfit_strong$linmod)[["meanalpha"]]
                                                        ),
                  from=0.0, to=1.0, col='tomato', add=TRUE, lwd=2 )


        # # Curve from approach 5
        # mycurve(  function(x) stress_exp_alpha( x,
        #                                         alpha=dplyr::select( data_tmp, meanalpha),
        #                                         coef(linearfit5$linmod_y0)[1], 
        #                                         coef(linearfit5$linmod_y0)[2], 
        #                                         coef(linearfit5$linmod_curve)[1], 
        #                                         coef(linearfit5$linmod_curve)[2] 
        #                                         ), 
        #           from=0.0, to=1.0, col='royalblue3', add=TRUE, lwd=2 )


        legend( "bottomright", c("approach I", "approach IV", "approach III"), lty=1, bty="n", lwd=2, col=c("springgreen3", "royalblue3", "tomato") )
  
      }
    }
    if (makepdf) dev.off()

  }

}