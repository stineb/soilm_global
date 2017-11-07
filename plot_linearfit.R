plot_linearfit <- function( linearfit, nlsfit=NULL ){

  par( las=1 )
  with( linearfit$df_flue0, plot( meanalpha, flue0, pch=16, xlab="AET/PET", ylab=expression(paste("fLUE"[0])), xlim=c(0,1) ) )
  abline( linearfit$linmod, col="black" )

  text( 0.3, 1.2, bquote( italic(R)^2 == .(format( summary( linearfit$linmod )$r.squared, digits = 2) ) ),  adj=0.0, cex=1 )
  cf <- coef(linearfit$linmod) %>% round( 2 )
  eq <- paste0( "y = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x " )
  text( 0.3, 1.1, eq, adj=0.0 )

  if (!is.null(nlsfit)){
    abline( a=coef(nlsfit)[[ "apar" ]], b=coef(nlsfit)[[ "bpar" ]], col="red" )
    eqfit <- paste0( "y = ", format( coef(nlsfit)[[ "apar" ]], digits=2 ), ifelse(sign(coef(nlsfit)[[ "bpar" ]])==1, " + ", " - "), format( abs(coef(nlsfit)[[ "bpar" ]]), digits=2 ), " x " )
    text( 0.3, 1.0, eqfit, adj=0.0, col="red" )
  }

}