library(dplyr)
library(ncdf4)
library(RColorBrewer)

##-----------------------------------------------
## Get all GPP time series
##-----------------------------------------------

## File names for TRENDY models
df_filnams <- read.csv("/Users/benjaminstocker/data/trendy/v5/trendy_s2_filnams_gpp.csv", stringsAsFactors = FALSE)

# modl           <- df_filnams$modl
# filnams        <- df_filnams$orig
modl           <- c()
filnams        <- c()
filnams_var    <- gsub( ".nc", "_sub_ann_var.nc", filnams )
filnams_relvar <- gsub( ".nc", "_sub_ann_relvar.nc", filnams )
path           <- paste0( "/Users/benjaminstocker/data/trendy/v5/", df_filnams$modl, "/S2/" )
filnams_var    <- filnams_var[ which(filnams!="") ]
filnams_relvar <- filnams_relvar[ which(filnams!="") ]
path           <- path[ which(filnams!="") ]
modl           <- modl[ which(filnams!="") ]
filpath_var    <- paste0( path, filnams_var )
filpath_relvar <- paste0( path, filnams_relvar )

filpath_var <- c( filpath_var, 
									"/Users/benjaminstocker/data/gpp_mte/MTE_var.nc", 
									"/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_var_s0_fapar3g_global.nc", 
									"/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_var_s1_fapar3g_global.nc"
									)
filpath_relvar <- c(  filpath_relvar, 
											"/Users/benjaminstocker/data/gpp_mte/MTE_relvar.nc", 
											"/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_relvar_s0_fapar3g_global.nc", 
											"/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_relvar_s1_fapar3g_global.nc"
											)
modl <- c( modl, "MTE", "Pmodel_S0", "Pmodel_S1")


## Load alpha (AET/PET) from P-model output to mask out values (some very high GPP interannual variance is not related to soil moisture)
nc <- nc_open( "~/data/pmodel_fortran_output/s0_fapar3g_global.a.alpha.nc" )
alpha <- ncvar_get( nc, varid="alpha" )
nc_close(nc)

## take mean alpha over the 30 years (1982-2011)
alpha <- apply( alpha, c(1,2), FUN = mean )

## define mask (considering only gridcells with non-zero water limitation)
mask <- ifelse( alpha<0.97, 1, NA )

relvar <- list()
var    <- list()
for (idx in seq(length(modl))){
# for (idx in 13:15){

	## Read GPP VARIANCE TRENDY file to complement dataframe
	if (file.exists(filpath_var[idx])){

		## read file
		nc <- nc_open( filpath_var[idx] )
		gpp <- try( ncvar_get( nc, varid="gpp" ) )
		nc_close(nc)

		## mask only dry-affected gridcells
		gpp <- gpp * mask

		if (class(gpp)!="try-error"){
		  ## transform values into a vector
		  vec <- c( gpp )
		  vec <- vec[ which( !is.na(vec) ) ]
		  if (length(vec)>0) var[[ modl[idx] ]] <- vec * 1e-3
		  # if (length(vec)>0) hist( var[[ modl[idx] ]], main=paste( modl[idx], "variance (kgC/m2/yr)"), breaks = 100 )
		}
	}

	## Read GPP RELATIVE VARIANCE TRENDY file to complement dataframe
	if (file.exists(filpath_relvar[idx])){

		## read file
		nc <- nc_open( filpath_relvar[idx] )
		gpp <- try( ncvar_get( nc, varid="gpp" ) )
		nc_close(nc)

		## mask only dry-affected gridcells
		gpp <- gpp * mask

		if (class(gpp)!="try-error"){
		  ## transform values into a vector
		  vec <- c( gpp )
		  vec <- vec[ which( !is.na(vec) ) ]
		  if (length(vec)>0) relvar[[ modl[idx] ]] <- vec
		  # if (length(vec)>0) hist( relvar[[ modl[idx] ]], main=paste( modl[idx], "relative variance (fraction)"), breaks = 100 )
		}
	}
	
}

save( relvar, file="data/relvar.Rdata" )
save( var, file="data/var.Rdata" )


##-----------------------------------------------
## Plot distribution of gridcell relative variance
##-----------------------------------------------
modls <- names(var)
modls_trendy <- modls[ -which( modls %in% c( "MTE", "Pmodel_S0", "Pmodel_S1" ) ) ]
cols_trendy <- c("springgreen3", "royalblue3", "tomato", "goldenrod", "orchid", "turquoise", "wheat", "darkslateblue", "grey70", "darkolivegreen" )[0:length(modls_trendy)]
cols <- c( cols_trendy, "red", "black", "black" ) 
lwds <- c(rep( 2, length(modls_trendy) ), rep(3,3) )
ltys <- c(rep( 1, length(modls_trendy) ), c(1,3,1) )
# cols_trendy <- brewer.pal( length(modls_trendy), "Set3" )
# cols_trendy <- colorRampPalette( brewer.pal(8, "Accent"))(length(modls_trendy))


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
par(las=1)
plot(  density( log10(var$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP absolute variance (kgC/m2/yr)")), ylim=c(0,1.2), xlim=c(-2,2), main="" )
lapply( seq(length(modls)), function(x) lines( density( log10(var[[ x ]]) ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topleft", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

# ## logarithmic x and y axis
# par(las=1)
# plot(  density( var$MTE ), col="red", type="n", lwd=2, xlab="GPP variance (kgC/m2/yr)", log="xy", xlim=c(0.5, 80), ylim=c(1e-4, 5e-1), main="" ) #
# lapply( seq(length(modls)), function(x) lines( density( var[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# # legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )


## Relative variance
##-----------------------------------------------
# ## straight axes
# par(las=1)
# plot(  density( relvar$MTE ), col="red", type="n", xlab="GPP relative variance (fraction)", ylim=c(0,0.25), xlim=c(0,60), main="" )
# lapply( seq(length(modls)), function(x) lines( density( relvar[[ x ]] ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
# legend( "topright", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

## logarithmic x axis (lognormal?)
par(las=1)
plot(  density( log10(relvar$MTE) ), col="red", type="n", xlab=expression(paste("log"[10]," of GPP relative variance (fraction)")), ylim=c(0,1.5), xlim=c(-2,2), main="" )
lapply( seq(length(modls)), function(x) lines( density( log10(relvar[[ x ]]) ), col=cols[x], lwd=lwds[x], lty=ltys[x] ) )
legend( "topleft", modls, col=cols, lwd=lwds, lty=ltys, bty = "n" )

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



