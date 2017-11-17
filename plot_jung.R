library(ncdf4)
library(dplyr)
source("~/.Rprofile")

overwrite <- FALSE
outfiln <- "data/ampl_jung.Rdata"

vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )

filpath_detr <- c(  paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_detr_s0_fapar3g_global.nc" ), 
                    paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_detr_s1_fapar3g_global.nc" )
                    )

filpath_nice <- c(  paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_nice_s0_fapar3g_global.nc" ), 
                    paste0( myhome, "/data/pmodel_fortran_output/pmodel_gpp_nice_s1_fapar3g_global.nc" )
                    )

modl <- c( "Pmodel_S0", "Pmodel_S1")

if (!file.exists(outfiln) || overwrite){

	ampl <- data.frame()
	ampl_u25 <- data.frame()

	for ( ires in 1:length(vec_res) ){

	  print(paste("resolution number", ires))
	  
		detr <- list()
		nice <- list()
		for (idx in seq(length(modl))){

		  ## Read detrended GPP for IAV
		  if (file.exists(filpath_detr[idx])){

		    ## read file
		    if (ires==1){
			    filn <- filpath_detr[idx]
			  } else {
			    filn <- gsub( "detr", paste0("detr_regr", sprintf("%02d", ires)), filpath_detr[idx] )
		    }
		    nc <- nc_open( filn )
		    detr[[ modl[idx] ]] <- ncvar_get( nc, varid="gpp" )
		    nc_close(nc)

		    ## read file
		    if (ires==1){
			    filn <- filpath_nice[idx]
			  } else {
			    filn <- gsub( "nice", paste0("nice_regr", sprintf("%02d", ires)), filpath_nice[idx] )
		    }
		    nc <- nc_open( filn )
		    nice[[ modl[idx] ]] <- ncvar_get( nc, varid="gpp" )
		    nc_close(nc)    

		  }
		  
		}

		## Get variance arrays for both simulations
		var <- list()
		relvar <- list()
		mean <- list()
		for (idx in seq(length(modl))){

			if (ires==length(vec_res)){

				var[[ modl[idx] ]]    <- var(  c( detr[[ modl[idx] ]] ), na.rm = TRUE )
				mean[[ modl[idx] ]]   <- mean( c( nice[[ modl[idx] ]] ), na.rm = TRUE )
				relvar[[ modl[idx] ]] <- var[[ modl[idx] ]] / mean[[ modl[idx] ]]

			} else {

				var[[ modl[idx] ]]    <- apply( detr[[ modl[idx] ]], c(1,2), FUN = var )
				mean[[ modl[idx] ]]   <- apply( nice[[ modl[idx] ]], c(1,2), FUN = mean )
				relvar[[ modl[idx] ]] <- var[[ modl[idx] ]] / mean[[ modl[idx] ]]

			}
		  
		}

		## Get amplification of relative variance, array and vector
		char_ires <- sprintf( "%02d", ires )
		ampl_arr <- relvar[[ "Pmodel_S1" ]] / relvar[[ "Pmodel_S0" ]]
		tmp_relvar <- c( ampl_arr )
		tmp_relvar <- tmp_relvar[ which( !is.na(tmp_relvar) ) ]

		## Get amplification of absolute variance, array and vector
		ampl_arr <- var[[ "Pmodel_S1" ]] / var[[ "Pmodel_S0" ]]
		tmp_var <- c( ampl_arr )
		tmp_var <- tmp_var[ which( !is.na(tmp_var) ) ]

		addrows <- data.frame( var=tmp_var, relvar=tmp_relvar, resnr=rep(ires, length(tmp_var)) )
		ampl <- rbind( ampl, addrows )	

		addrows <- addrows %>% filter( relvar > quantile( relvar, probs=0.75 ) )
		ampl_u25 <- rbind( ampl_u25, addrows )	

	}

	save( ampl, file=outfiln )

	# boxplot( var ~ resnr, data=ampl )
	# boxplot( relvar ~ resnr, data=ampl_u25 )


	save( ampl_agg, file="data/ampl_agg_jung.Rdata" )

} else {

	load( outfiln )

}

## aggregate by resolution
ampl_agg <- ampl %>%  group_by( resnr ) %>% 
												summarise(  var_mean = mean( var, na.rm=TRUE ), 
																		var_median = median( var, na.rm=TRUE ),
																		var_q01 = quantile( var, probs=0.01 ),
																		var_q05 = quantile( var, probs=0.05 ),
																		var_q10 = quantile( var, probs=0.10 ),
																		var_q90 = quantile( var, probs=0.90 ),
																		var_q95 = quantile( var, probs=0.95 ),
																		var_q99 = quantile( var, probs=0.99 ),
																		var_q25 = quantile( var, probs=0.25 ),
																		var_q75 = quantile( var, probs=0.75 ),

																		relvar_mean = mean( relvar, na.rm=TRUE ), 
																		relvar_median = median( relvar, na.rm=TRUE ),
																		relvar_q01 = quantile( relvar, probs=0.01 ),
																		relvar_q05 = quantile( relvar, probs=0.05 ),
																		relvar_q10 = quantile( relvar, probs=0.10 ),
																		relvar_q90 = quantile( relvar, probs=0.90 ),
																		relvar_q95 = quantile( relvar, probs=0.95 ),
																		relvar_q99 = quantile( relvar, probs=0.99 ),
																		relvar_q25 = quantile( relvar, probs=0.25 ),
																		relvar_q75 = quantile( relvar, probs=0.75 )
																		)

## Plot scale dependence of soil moisture effect on GPP interannual variance
par(las=1)
with( ampl_agg, plot( resnr, relvar_mean, type="l", col="black", lwd=2, ylim=c(0,6.5), xaxt = "n", xlab="spatial resolution (degrees)", ylab="ampl. of rel. var. of ann. GPP" ) )
axis( 1, at=seq(length(vec_res)), labels=as.character( vec_res ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q01, relvar_q99), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q05, relvar_q95), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q10, relvar_q90), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q25, relvar_q75), col=add_alpha("tomato2", 0.25), border = NA ) )
abline( h = 1.0, lty=3 )





