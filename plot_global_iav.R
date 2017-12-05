library(dplyr)
library(tibble)
library(ncdf4)

##-----------------------------------------------
## Get all GPP time series
##-----------------------------------------------

# ## File names for TRENDY models
# filnams <- read.csv( paste0( myhome, "/data/trendy/v5/trendy_s2_filnams_gpp.csv" ) )

# # rm("df_gpp")

# i <- 0
# for (modl in filnams$modl){

# 	## Read TRENDY file to complement dataframe
#   i <- i + 1
# 	df_tmp <- try( read.csv( paste0( "data/", modl, "_globaltotal.csv" ) ) %>% setNames( c( "year", modl ) ) )
	
# 	if (class(df_tmp)!="try-error"){
# 		if (i==1){
# 			df_gpp <- df_tmp
# 		} else {
# 			df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )
# 		}		
# 	}

# }

## MTE
## read file
nc <- nc_open( paste0( myhome, "/data/gpp_mte/gpp_mte_GLOB30yr.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
time <- 1982:2011
nc_close(nc)
df_gpp <- tibble( year=time, MTE=gpp )
# df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )

## MTE FLUXCOM
nc <- nc_open( paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_GLOB30yr.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
time <- 1982:2011
nc_close(nc)
df_tmp <- tibble( year=time, MTE_FLUXCOM=gpp )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )


## P-model, s0 simulation (without soil moisture stress)
nc <- nc_open( paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_GLOB30yr.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
time <- 1982:2011
nc_close(nc)
df_tmp <- tibble( year=time, pmodel_s0=gpp )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )

## P-model, s1 simulation (with soil moisture stress)
nc <- nc_open( paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1_GLOB30yr.nc" ) )
gpp <- ncvar_get( nc, varid="gpp" )
time <- 1982:2011
nc_close(nc)
df_tmp <- tibble( year=time, pmodel_s1=gpp )
df_gpp <- df_gpp %>% left_join( df_tmp, by="year" )

## Drop columns if all values are NA
df_gpp <- df_gpp %>% dplyr::select( which( colMeans(is.na(.)) < 0.1 ) )

## Rename because dplyr can't handle '-'
df_gpp <- df_gpp %>% rename( "LPJ_GUESS"="LPJ-GUESS", "LPX_Bern"="LPX-Bern" )

## Get mean GPP across all years for each model 
df_gpp_mean <- df_gpp %>% summarise(  
                                      # CABLE=mean(CABLE), 
																			# CLM=mean(CLM), 
																			# ISAM=mean(ISAM), 
																			# JSBACH=mean(JSBACH), 
																			# LPJ_GUESS=mean(LPJ_GUESS), 
																			# LPX_Bern=mean(LPX_Bern), 
																			# ORCHIDEE=mean(ORCHIDEE), 
																			# SDGVM=mean(SDGVM), 
																			# VEGAS=mean(VEGAS), 
																			# VISIT=mean(VISIT), 
																			MTE=mean(MTE), 
																			MTE_FLUXCOM=mean(MTE_FLUXCOM), 
																			pmodel_s0=mean(pmodel_s0), 
																			pmodel_s1=mean(pmodel_s1) 
																			)

df_gpp_var <- df_gpp %>% summarise( 
                                    # CABLE=var(CABLE), 
																		# CLM=var(CLM), 
																		# ISAM=var(ISAM), 
																		# JSBACH=var(JSBACH), 
																		# LPJ_GUESS=var(LPJ_GUESS), 
																		# LPX_Bern=var(LPX_Bern), 
																		# ORCHIDEE=var(ORCHIDEE), 
																		# SDGVM=var(SDGVM), 
																		# VEGAS=var(VEGAS), 
																		# VISIT=var(VISIT), 
																		MTE=var(MTE), 
																		MTE_FLUXCOM=var(MTE_FLUXCOM), 
																		pmodel_s0=var(pmodel_s0), 
																		pmodel_s1=var(pmodel_s1) 
																		)


# ###----------------------------------------------
# ## Get linear model for each time series (model) and extract residuals (effectively detrending the time series)
# ##-----------------------------------------------
# modls <- names( df_gpp )[ -which( names(df_gpp)=="year" ) ]
# df_detr <- df_gpp
# df_detr_rel <- df_gpp
# 
# for (modl in modls){
# 	df_tmp <- df_gpp %>% select( year, modl )
# 	lm_tmp <- lm( df_tmp[[modl]] ~ df_tmp$year )
# 	df_detr[[ modl ]] <- lm_tmp$residuals
# 	df_detr_rel[[ modl ]] <- 100 * lm_tmp$residuals / df_gpp_mean[[ modl ]]
# }

##-----------------------------------------------
## Plot 
##-----------------------------------------------
modls_trendy <- modls[ -which( modls %in% c( "MTE", "MTE_FLUXCOM", "pmodel_s0", "pmodel_s1" ) ) ]
cols_trendy <- c("springgreen3", "royalblue3", "tomato", "goldenrod", "orchid", "turquoise", "wheat", "darkslateblue", "grey70", "darkolivegreen" ) 
cols <- c( cols_trendy, "red", "red", "black", "black" ) 
lwds <- c(rep( 2, length(modls_trendy) ), rep(3,4) )
ltys <- c(rep( 1, length(modls_trendy) ), c(1,3,1,3) )

##-----------------------------------------------
## Plot IAV: distribution of annual global GPP anomalies from mean trend
##-----------------------------------------------
par(las=1)
plot( density( df_detr$MTE ), xlim=c(-5,5), ylim=c(0,0.5), main="Interannual variability (1982-2011)", xlab=expression(paste("global total GPP anomaly (PgC yr"^{-1}, ")")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr[[ x + 1 ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

##-----------------------------------------------
## Plot relative IAV: distribution of annual global GPP anomalies (normalised by mean GPP) from mean trend
##-----------------------------------------------
pdf( "fig/relative_IAV.pdf")
par(las=1)
plot( density( df_detr_rel$MTE ), xlim=c(-5,5), ylim=c(0,0.75), main="Interannual variability, relative to mean (1982-2011)", xlab=expression(paste("global total GPP anomaly (%)")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr_rel[[ x + 1 ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )
dev.off()




##----------------------------------------------
## Short time series (10 yr)
##-----------------------------------------------
## P-model, s0 simulation (without soil moisture stress)
df_gpp_10y <- try( read.csv( paste0( myhome, "data/pmodel_fortran_output/pmodel_gpp_globaltotals0_fapar3g_global_10y.csv" ) ) %>% setNames( c( "year", "pmodel_s0" ) ) )

## P-model, s1 simulation (with soil moisture stress)
df_tmp <- try( read.csv( paste0( myhome, "data/pmodel_fortran_output/pmodel_gpp_globaltotals1_fapar3g_global_10y.csv" ) ) %>% setNames( c( "year", "pmodel_s1" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )

## MTE
df_tmp <- try( read.csv( paste0( myhome, "/data/gpp_mte/gpp_mte_globaltotal_10y.csv" ) ) %>% setNames( c( "year", "MTE" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )

## MTE FLUXCOM
df_tmp <- try( read.csv( paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_globaltotal_10y.csv" ) ) %>% setNames( c( "year", "MTE_FLUXCOM" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )

## VPM
df_tmp <- try( read.csv( paste0( myhome, "/data/gpp_vpm/gpp_vpm_globaltotal_10y.csv" ) ) %>% setNames( c( "year", "VPM" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )

## BESS
df_tmp <- try( read.csv( paste0( myhome, "/data/gpp_bess/gpp_bess_globaltotal_10y.csv" ) ) %>% setNames( c( "year", "BESS" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )

## MODIS
df_tmp <- try( read.csv( paste0( myhome, "/data/gpp_modis/gpp_modis_globaltotal_10y.csv" ) ) %>% setNames( c( "year", "MODIS" ) ) )
df_gpp_10y <- df_gpp_10y %>% left_join( df_tmp, by="year" )


## Get mean GPP across all years for each model 
df_gpp_mean_10y <- df_gpp_10y %>% summarise(  MTE=mean(MTE), 
																							MTE_FLUXCOM=mean(MTE_FLUXCOM), 
																							VPM=mean(VPM), 
																							BESS=mean(BESS), 
																							MODIS=mean(MODIS), 
																							pmodel_s0=mean(pmodel_s0), 
																							pmodel_s1=mean(pmodel_s1) 
																							)

df_gpp_var_10y <- df_gpp_10y %>% summarise( MTE=var(MTE), 
                                            MTE_FLUXCOM=var(MTE_FLUXCOM), 
                                            VPM=var(VPM), 
                                            BESS=var(BESS), 
                                            MODIS=var(MODIS), 
                                            pmodel_s0=var(pmodel_s0), 
                                            pmodel_s1=var(pmodel_s1) 
)

###----------------------------------------------
## Get linear model for each time series (model) and extract residuals (effectively detrending the time series)
##-----------------------------------------------
modls <- names( df_gpp_10y )[ -which( names(df_gpp_10y)=="year" ) ]
df_detr_10y <- df_gpp_10y
df_detr_rel_10y <- df_gpp_10y

for (modl in modls){
	df_tmp <- df_gpp_10y %>% select( year, modl )
	lm_tmp <- lm( df_tmp[[modl]] ~ df_tmp$year )
	df_detr_10y[[ modl ]] <- lm_tmp$residuals
	df_detr_rel_10y[[ modl ]] <- 100 * lm_tmp$residuals / df_gpp_mean_10y[[ modl ]]
}

##-----------------------------------------------
## Plot 
##-----------------------------------------------
cols <- c("springgreen3", "royalblue3", "tomato", "goldenrod", "orchid", "turquoise", "wheat" ) 
lwds <- rep( 2, length(modls) )
ltys <- rep( 1, length(modls) )

##-----------------------------------------------
## Plot IAV: distribution of annual global GPP anomalies from mean trend
##-----------------------------------------------
par(las=1)
plot( density( df_detr_10y$MTE_FLUXCOM ), xlim=c(-5,5), ylim=c(0,2), main="Interannual variability (2001-2011)", xlab=expression(paste("global total GPP anomaly (PgC yr"^{-1}, ")")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr_10y[[ x + 1 ]], adjust=2 ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

##-----------------------------------------------
## Plot relative IAV: distribution of annual global GPP anomalies (normalised by mean GPP) from mean trend
##-----------------------------------------------
pdf( "fig/relative_IAV_10y.pdf")
par(las=1)
plot( density( df_detr_rel_10y$MTE_FLUXCOM ), xlim=c(-5,5), ylim=c(0,2), main="Interannual variability, relative to mean (2001-2011)", xlab=expression(paste("global total GPP anomaly (%)")), type="n" )
lapply( seq(length(modls)), function(x) lines( density( df_detr_rel_10y[[ x + 1 ]], adjust=2 ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )
dev.off()




