plot_fit_fvar_vs_time <- function( linearfit1, linearfit2, linearfit3, linearfit5, ddf=NULL, nice_agg=NULL , makepdf=FALSE){

  require(dplyr)
  require(lubridate)

  if (!is.null(ddf)){

    ##-----------------------------------------------
    ## Time series: fLUE and fLUEest from approaches I, II, and III
    ##-----------------------------------------------
    filn <- "fig/flue_est_per_site.pdf" 
    if (makepdf) print( paste( "plotting fLUE and fLUEest vs. time for each site into file ", filn, "..." ) )
    if (makepdf) pdf( filn, width = 10, height = 6 )
    for (sitename in linearfit2$data$mysitename){

      df_tmp <- dplyr::filter(nice_agg, mysitename==sitename)

      if (nrow(df_tmp)>0){
        
        par(las=1)
        plot(  df_tmp$date, df_tmp[[ "fvar" ]], type="l", xlab="time", ylab="fLUE", col="black", ylim = c(0.2, 1.2) )
        lines( df_tmp$date, df_tmp$flue_est_1, col="springgreen3" )
        lines( df_tmp$date, df_tmp$flue_est_5, col="royalblue3" )
        lines( df_tmp$date, df_tmp$flue_est_3, col="tomato" )
        title( sitename )
        legend( "bottomright", c("approach I", "approach II", "approach III"), lty=1, bty="n", lwd=2, col=c("springgreen3", "royalblue3", "tomato") )
      }
    }
    if (makepdf) dev.off()

  }

}