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
cols <- c( cols_trendy, "red", "black", "black" ) 
lwds <- c(rep( 2, length(modls_trendy) ), rep(3,3) )
ltys <- c(rep( 1, length(modls_trendy) ), rep(1,2), 3 )

par(las=1)
plot( density( df_detr$CLM ), xlim=c(-5,5), ylim=c(0,0.5), main="Interannual variability (1982-2011)", xlab=expression(paste("global total GPP anomaly (PgC yr"^{-1}, ")")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr[[ x + 1 ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )






