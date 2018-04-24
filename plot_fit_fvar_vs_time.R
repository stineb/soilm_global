plot_fit_fvar_vs_time <- function( linearfit1, linearfit_mid, linearfit_strong, ddf=NULL, nice_agg=NULL , makepdf=FALSE){

  require(dplyr)
  require(lubridate)
  require(zoo)

  if (!is.null(ddf)){
    ##-----------------------------------------------
    ## Time series: fLUE and fLUEest from approaches I, II, and III
    ##-----------------------------------------------
    filn <- "fig/flue_est_per_site.pdf" 
    if (makepdf) print( paste( "plotting fLUE and fLUEest vs. time for each site into file ", filn, "..." ) )
    if (makepdf) pdf( filn, width = 10, height = 6 )
    for (sitename in linearfit_mid$data$mysitename){

      df_tmp <- dplyr::filter(nice_agg, mysitename==sitename)

      if (nrow(df_tmp)>0){
        
        par(las=1)
        plot(  df_tmp$date, df_tmp[[ "fvar" ]], type="l", xlab="time", ylab="fLUE", col="black", ylim = c(0.0, 1.2) )
        # lines( df_tmp$date, df_tmp[[ "fvar_smooth" ]], col="black" )
        lines( df_tmp$date, df_tmp$flue_est_I, col="springgreen3" )
        lines( df_tmp$date, df_tmp$flue_est_IV, col="royalblue3" )
        lines( df_tmp$date, df_tmp$flue_est_II, col="tomato" )
        title( sitename )
        legend( "bottomright", c("approach I", "approach IV", "approach III"), lty=1, bty="n", lwd=2, col=c("springgreen3", "royalblue3", "tomato") )
        abline( h=0.2, lty=3 )
        
      }
    }
    if (makepdf) dev.off()

  }

}