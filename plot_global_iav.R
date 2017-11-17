library(dplyr)

##-----------------------------------------------
## Get all GPP time series
##-----------------------------------------------

## File names for TRENDY models
filnams <- read.csv( paste0( myhome, "/data/trendy/v5/trendy_s2_filnams_gpp.csv" ) )

# rm("df_gpp")

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

## Rename because dplyr can't handle '-'
df_gpp <- df_gpp %>% rename( "LPJ_GUESS"="LPJ-GUESS", "LPX_Bern"="LPX-Bern" )

## Get mean GPP across all years for each model 
df_gpp_mean <- df_gpp %>% summarise(  CABLE=mean(CABLE), 
																			CLM=mean(CLM), 
																			ISAM=mean(ISAM), 
																			JSBACH=mean(JSBACH), 
																			LPJ_GUESS=mean(LPJ_GUESS), 
																			LPX_Bern=mean(LPX_Bern), 
																			ORCHIDEE=mean(ORCHIDEE), 
																			SDGVM=mean(SDGVM), 
																			VEGAS=mean(VEGAS), 
																			VISIT=mean(VISIT), 
																			MTE=mean(MTE), 
																			pmodel_s0=mean(pmodel_s0), 
																			pmodel_s1=mean(pmodel_s1) 
																			)

###----------------------------------------------
## Get linear model for each time series (model) and extract residuals (effectively detrending the time series)
##-----------------------------------------------
modls <- names( df_gpp )[ -which( names(df_gpp)=="year" ) ]
df_detr <- df_gpp
df_detr_rel <- df_gpp

for (modl in modls){
	df_tmp <- df_gpp %>% select( year, modl )
	lm_tmp <- lm( df_tmp[[modl]] ~ df_tmp$year )
	df_detr[[ modl ]] <- lm_tmp$residuals
	df_detr_rel[[ modl ]] <- 100 * lm_tmp$residuals / df_gpp_mean[[ modl ]]
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

par(las=1)
plot( density( df_detr_rel$CLM ), xlim=c(-5,5), ylim=c(0,0.75), main="Interannual variability, relative to mean (1982-2011)", xlab=expression(paste("global total GPP anomaly (%)")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr_rel[[ x + 1 ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )






