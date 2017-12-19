library(dplyr)
library(ncdf4)
library(RColorBrewer)
source("~/.Rprofile")

##-----------------------------------------------
## Get all GPP time series
##-----------------------------------------------

## File names for TRENDY models
df_filnams <- read.csv( paste0( myhome, "/data/trendy/v5/trendy_s2_filnams_gpp.csv"), stringsAsFactors = FALSE)

modls_trendy    <- df_filnams$modl
filnams         <- df_filnams$orig
# modls_trendy       <- c()
# filnams        <- c()
filnams_var_trendy    <- gsub( ".nc", "_VAR.nc", filnams )
filnams_relvar_trendy <- gsub( ".nc", "_RELVAR.nc", filnams )
path_trendy           <- paste0(  myhome, "/data/trendy/v5/", modls_trendy, "/S2/" )
filnams_var_trendy    <- filnams_var_trendy[ which(filnams!="") ]
filnams_relvar_trendy <- filnams_relvar_trendy[ which(filnams!="") ]
path_trendy           <- path_trendy[ which(filnams!="") ]
modls_trendy          <- modls_trendy[ which(filnams!="") ]
filpath_var_trendy    <- paste0( path_trendy, filnams_var_trendy )
filpath_relvar_trendy <- paste0( path_trendy, filnams_relvar_trendy )

## RS-models with 30 year data availability
modls_30y <- c( "MTE", "MTE FLUXCOM", "P-model s0", "P-model s1")
filpath_var_30y <- c(  paste0( myhome, "/data/gpp_mte/gpp_mte_VAR.nc"), 
											 paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_VAR.nc"), 
											 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_VAR.nc"), 
											 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1_VAR.nc")
											)
filpath_relvar_30y <- c(   paste0( myhome, "/data/gpp_mte/gpp_mte_RELVAR.nc"), 
													 paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_RELVAR.nc"), 
													 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_RELVAR.nc"), 
													 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1_RELVAR.nc")
													)

## RS-models with 10 year data availability
modls_10y <- c( "MTE", "MTE FLUXCOM", "P-model s0", "P-model s1", "MOD17AH2", "BESS", "VPM")
filpath_var_10y <- c(  paste0( myhome, "/data/gpp_mte/gpp_mte_VAR20XX.nc"), 
											 paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_VAR20XX.nc"), 
											 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_VAR20XX.nc"), 
											 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1_VAR20XX.nc"),
											 paste0( myhome, "/data/gpp_modis/gpp_modis_VAR20XX.nc"),
											 paste0( myhome, "/data/gpp_bess/gpp_bess_VAR20XX.nc"),
											 paste0( myhome, "/data/gpp_vpm/gpp_vpm_VAR20XX.nc")
											)
filpath_relvar_10y <- c(   paste0( myhome, "/data/gpp_mte/gpp_mte_RELVAR20XX.nc"), 
													 paste0( myhome, "/data/gpp_mte/gpp_mte_fluxcom_RELVAR20XX.nc"), 
													 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_RELVAR20XX.nc"), 
													 paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1_RELVAR20XX.nc"),
													 paste0( myhome, "/data/gpp_modis/gpp_modis_RELVAR20XX.nc"),
													 paste0( myhome, "/data/gpp_bess/gpp_bess_RELVAR20XX.nc"),
													 paste0( myhome, "/data/gpp_vpm/gpp_vpm_RELVAR20XX.nc")
													)

modl <- list( len_30y=modls_30y, len_10y=modls_10y, trendy=modls_trendy )
filpath_var <- list( len_30y=filpath_var_30y, len_10y=filpath_var_10y, trendy=filpath_var_trendy )
filpath_relvar <- list( len_30y=filpath_relvar_30y, len_10y=filpath_relvar_10y, trendy=filpath_relvar_trendy )

## Load alpha (AET/PET) from P-model output to mask out values (some very high GPP interannual variance is not related to soil moisture)
nc <- nc_open(  paste0( myhome, "/data/pmodel_fortran_output/s0_fapar3g_global.a.alpha.nc") )
alpha <- ncvar_get( nc, varid="alpha" )
nc_close(nc)

## take mean alpha over the 30 years (1982-2011)
alpha <- apply( alpha, c(1,2), FUN = mean )

## define mask (considering only gridcells with non-zero water limitation)
mask <- ifelse( alpha<0.97, 1, NA )

relvar <- list()
var    <- list()

for (iscal in c("len_30y", "len_10y", "trendy")){

	for (idx in seq( length( modl[[ iscal ]] ) )){

		## Read GPP VARIANCE TRENDY file to complement dataframe
		if (file.exists(filpath_var[[ iscal ]][idx])){

			## read file
			nc <- nc_open( filpath_var[[ iscal ]][idx] )
			gpp <- try( ncvar_get( nc, varid="gpp" ) )
			nc_close(nc)

			## mask only dry-affected gridcells
			if (iscal!="trendy") gpp <- gpp * mask

			if (class(gpp)!="try-error"){
			  ## transform values into a vector
			  vec <- c( gpp )
			  vec <- vec[ which( !is.na(vec) ) ]
			  if (length(vec)>0) var[[ iscal ]][[ modl[[ iscal ]][idx] ]] <- vec * 1e-3
			  # if (length(vec)>0) hist( var[[ modl[[ iscal ]][idx] ]], main=paste( modl[[ iscal ]][idx], "variance (kgC/m2/yr)"), breaks = 100 )
			}
		} else {
			print( paste( "file does not exist:", filpath_var[[ iscal ]][idx] ))
		}

		## Read GPP RELATIVE VARIANCE TRENDY file to complement dataframe
		if (file.exists(filpath_relvar[[ iscal ]][idx])){

			## read file
			nc <- nc_open( filpath_relvar[[ iscal ]][idx] )
			gpp <- try( ncvar_get( nc, varid="gpp" ) )
			nc_close(nc)

			## mask only dry-affected gridcells
			if (iscal!="trendy") gpp <- gpp * mask

			if (class(gpp)!="try-error"){
			  ## transform values into a vector
			  vec <- c( gpp )
			  vec <- vec[ which( !is.na(vec) ) ]
			  if (length(vec)>0) relvar[[ iscal ]][[ modl[[ iscal ]][idx] ]] <- vec
			  # if (length(vec)>0) hist( relvar[[ modl[[ iscal ]][idx] ]], main=paste( modl[[ iscal ]][idx], "relative variance (fraction)"), breaks = 100 )
			}
		} else {
			print( paste( "file does not exist:", filpath_relvar[[ iscal ]][idx] ))
		}
		
	}

}

save( relvar, file="data/relvar.Rdata" )
save( var, file="data/var.Rdata" )


##-----------------------------------------------
## Plot distribution of gridcell relative variance
##-----------------------------------------------
cols <- list( len_30y=c( "red", "orchid", "black", "black" ), 
							len_10y=c( "red", "orchid", "black", "black", "springgreen3", "royalblue3", "goldenrod" ),
							trendy=c( "springgreen3", "royalblue3", "tomato", "goldenrod", "orchid", "turquoise", "wheat", "darkslateblue", "grey70", "darkolivegreen" )[1:length(var$trendy)] 
							)
lwds <- list( len_30y=rep(2,4), len_10y=rep(2,7), trendy=rep(2,length(var$trendy)) )
ltys <- list( len_30y=c(1,1,3,1), len_10y=c(1,1,3,1,1,1,1), trendy=rep(1,length(var$trendy)) )

##-----------------------------------------------
## Absolute variance
##-----------------------------------------------
# ## straight axes
# par(las=1)
# plot(  density( var$MTE ), col="red", type="n", xlab="GPP variance (kgC/m2/yr)", xlim=c(0,40), ylim=c(0,0.4), main="" )
# lapply( seq(length(modls)), function(x) lines( density( var[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

# ## logarithmic y axis
# par(las=1)
# plot(  density( var$MTE ), col="red", type="n", lwd=2, xlab="GPP variance (kgC/m2/yr)", log="y", xlim=c(0,40), ylim=c(1e-4, 5e-1) )
# lapply( seq(length(modls)), function(x) lines( density( var[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

## logarithmic x axis (lognormal?)
## RS-models, 30y
par(las=1)
plot(  density( log10(var$len_30y$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP absolute variance (kgC/m2/yr)")), ylim=c(0,1.2), xlim=c(-2,2.5), main="variance, RS, 1982-2011" )
lapply( seq(length(modl$len_30y)), function(x) lines( density( log10(var$len_30y[[ x ]]) ), col=cols$len_30y[x], lwd=lwds$len_30y[x], lty=ltys$len_30y[x] ) )
legend( "topleft", modl$len_30y, col=cols$len_30y, lwd=lwds$len_30y, lty=ltys$len_30y, bty = "n" )

## RS-models, 10y
plot(  density( log10(var$len_10y$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP absolute variance (kgC/m2/yr)")), ylim=c(0,1.2), xlim=c(-2,2.5), main="variance, RS, 2001-2011" )
lapply( seq(length(modl$len_10y)), function(x) lines( density( log10(var$len_10y[[ x ]]) ), col=cols$len_10y[x], lwd=lwds$len_10y[x], lty=ltys$len_10y[x] ) )
legend( "topleft", modl$len_10y, col=cols$len_10y, lwd=lwds$len_10y, lty=ltys$len_10y, bty = "n" )

## TRENDY models, 30y
plot(  density( log10(var$trendy$CLM) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP relative variance (fraction)")), ylim=c(0,1.2), xlim=c(-2,2.5), main="variance, TRENDY, 1982-2011" )
lapply( seq(length(var$trendy)), function(x) lines( density( log10( var$trendy[[ x ]]) ), col=cols$trendy[x], lwd=lwds$trendy[x], lty=ltys$trendy[x] ) )
legend( "topleft", modl$trendy, col=cols$trendy, lwd=lwds$trendy, lty=ltys$trendy, bty = "n" )
## add P-model
lapply( 3:4, function(x) lines( density( log10(var$len_30y[[ x ]]) ), col=cols$len_30y[x], lwd=lwds$len_30y[x], lty=ltys$len_30y[x] ) )
legend( "topright", modl$len_30y[3:4], col=cols$len_30y[3:4], lwd=lwds$len_30y[3:4], lty=ltys$len_30y[3:4], bty = "n" )


# ## logarithmic x and y axis
# par(las=1)
# plot(  density( var$MTE ), col="red", type="n", lwd=2, xlab="GPP variance (kgC/m2/yr)", log="xy", xlim=c(0.5, 80), ylim=c(1e-4, 5e-1), main="" ) #
# lapply( seq(length(modls)), function(x) lines( density( var[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# # legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )


##-----------------------------------------------
## Relative variance
##-----------------------------------------------
# ## straight axes
# par(las=1)
# plot(  density( relvar$MTE ), col="red", type="n", xlab="GPP relative variance (fraction)", ylim=c(0,0.25), xlim=c(0,60), main="" )
# lapply( seq(length(modls)), function(x) lines( density( relvar[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

## logarithmic x axis (lognormal?)
## RS-models, 30y
plot(  density( log10(relvar$len_30y$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP relative variance (fraction)")), ylim=c(0,1.5), xlim=c(-2,3), main="relative variance, RS, 1982-2011" )
lapply( seq(length(modl$len_30y)), function(x) lines( density( log10( relvar$len_30y[[ x ]]) ), col=cols$len_30y[x], lwd=lwds$len_30y[x], lty=ltys$len_30y[x] ) )
legend( "topleft", modl$len_30y, col=cols$len_30y, lwd=lwds$len_30y, lty=ltys$len_30y, bty = "n" )

## RS-models, 10y
plot(  density( log10(relvar$len_10y$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP relative variance (fraction)")), ylim=c(0,1.5), xlim=c(-2,3), main="relative variance, RS, 2001-2011" )
lapply( seq(length(modl$len_10y)), function(x) lines( density( log10( relvar$len_10y[[ x ]]) ), col=cols$len_10y[x], lwd=lwds$len_10y[x], lty=ltys$len_10y[x] ) )
legend( "topleft", modl$len_10y, col=cols$len_10y, lwd=lwds$len_10y, lty=ltys$len_10y, bty = "n" )

## TRENDY models, 30y
plot(  density( log10(relvar$trendy$CLM) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP relative variance (fraction)")), ylim=c(0,1.5), xlim=c(-2,3), main="relative variance, TRENDY, 1982-2011" )
lapply( seq(length(var$trendy)), function(x) lines( density( log10( relvar$trendy[[ x ]]) ), col=cols$trendy[x], lwd=lwds$trendy[x], lty=ltys$trendy[x] ) )
legend( "topleft", modl$trendy, col=cols$trendy, lwd=lwds$trendy, lty=ltys$trendy, bty = "n" )
## add P-model
lapply( 3:4, function(x) lines( density( log10( relvar$len_30y[[ x ]]) ), col=cols$len_30y[x], lwd=lwds$len_30y[x], lty=ltys$len_30y[x] ) )
legend( "topright", modl$len_30y[3:4], col=cols$len_30y[3:4], lwd=lwds$len_30y[3:4], lty=ltys$len_30y[3:4], bty = "n" )

# par(las=1)
# plot(  density( (relvar$MTE) ), col="red", type="n", xlab="GPP relative variance (fraction)", main="", log="x" )
# lapply( seq(length(modls)), function(x) lines( density( (relvar[[ x ]]) ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

# ## logarithmic y axis
# par(las=1)
# plot(  density( relvar$MTE ), col="red", type="n", lwd=2, xlab="GPP relative variance (fraction)", log="y", ylim=c(1e-6, 5e-1), xlim=c(0,60), main="" )
# lapply( seq(length(modls)), function(x) lines( density( relvar[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

# ## logarithmic x and y axis
# par(las=1)
# plot(  density( relvar$MTE ), col="red", type="n", lwd=2, xlab="GPP relative variance (fraction)", log="xy", ylim=c(1e-4, 2e-1), xlim=c(0.5,60), main="" ) #,
# lapply( seq(length(modls)), function(x) lines( density( relvar[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# # legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )



