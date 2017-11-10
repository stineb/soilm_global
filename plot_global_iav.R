library(dplyr)

##-----------------------------------------------
## Get all GPP time series
##-----------------------------------------------

## File names for TRENDY models
filnams <- read.csv("/Users/benjaminstocker/data/trendy/v5/trendy_s2_filnams_gpp.csv")

rm("df_gpp")

i <- 0
for (modl in filnams$modl){

	## Read TRENDY file to complement dataframe
  i <- i + 1
	df_tmp <- try( read.csv( paste0( "data/", modl, "_globaltotal.csv" ) ) %>% setNames( c( "year", modl ) ) )
	
	if (class(df_tmp)!="try-error"){
		if (i==1){
			df_gpp <- df_tmp
		} else {
			df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )
		}		
	}

}

## MTE
df_tmp <- try( read.csv( paste0( "data/mte_globaltotal.csv" ) ) %>% setNames( c( "year", "MTE" ) ) )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )

## P-model, s0 simulation (without soil moisture stress)
df_tmp <- try( read.csv( paste0( "data/pmodel_gpp_globaltotals0_fapar3g_global.csv" ) ) %>% setNames( c( "year", "pmodel_s0" ) ) )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )

## P-model, s1 simulation (with soil moisture stress)
df_tmp <- try( read.csv( paste0( "data/pmodel_gpp_globaltotals1_fapar3g_global.csv" ) ) %>% setNames( c( "year", "pmodel_s1" ) ) )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )


## Drop columns if all values are NA
df_gpp <- df_gpp %>% select( which( colMeans(is.na(.)) < 0.1 ) )


##-----------------------------------------------
## Get linear model for each time series (model) and extract residuals (effectively detrending the time series)
##-----------------------------------------------
modls <- names( df_gpp )[ -which( names(df_gpp)=="year" ) ]
df_detr <- df_gpp

for (modl in modls){
	df_tmp <- df_gpp %>% select( year, modl )
	lm_tmp <- lm( df_tmp[[modl]] ~ df_tmp$year )
	df_detr[[ modl ]] <- lm_tmp$residuals
}

##-----------------------------------------------
## Plot distribution of annual global GPP anomalies from mean trend
##-----------------------------------------------
modls_trendy <- modls[ -which( modls %in% c( "MTE", "pmodel_s0", "pmodel_s1" ) ) ]
cols_trendy <- c("springgreen3", "royalblue3", "tomato", "goldenrod", "orchid", "turquoise", "wheat", "darkslateblue", "grey70", "darkolivegreen" ) 
df_col_trendy <- data.frame( modl=modls_trendy, col=cols_trendy )

par(las=1)
plot( density( df_detr$CABLE ), xlim=c(-5,5), ylim=c(0,0.5), xlab="global total GPP anomaly (PgC/yr)", main="", col=as.character( filter( df_col_trendy, modl=="CABLE" )$col ), lwd=2 )
for (modl in modls_trendy){
	print( paste( "model:",modl, "col:", as.character( df_col_trendy[ which( df_col_trendy$modl==modl ), 2 ] ) ) ) 
	lines( density( df_detr[[ modl ]] ), col=as.character( df_col_trendy[ which( df_col_trendy$modl==modl ), 2 ] ), lwd=2 )
}
legend( "topleft", as.character(df_col_trendy$modl), col=as.character(df_col_trendy$col), bty = "n", lwd = 2, lty = 1 )

lines( density( df_detr$MTE ), col="red", lwd=3 )
legend( "topleft", "MTE", col="red", bty = "n", lwd = 3, lty = 1, inset = c(0,0.35) )

lines( density( df_detr$pmodel_s0 ), col="black", lwd=3 )
lines( density( df_detr$pmodel_s1 ), col="black", lwd=3, lty = 2 )
legend( "topleft", c("P-model, unstressed", "P-model, stressed"), col="black", bty = "n", lwd = 3, lty = c(1,2), inset = c(0,0.4) )







