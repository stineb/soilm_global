library(ncdf4)
library(dplyr)
library(abind)
source("~/.Rprofile")

overwrite <- TRUE
outfiln <- "data/ampl_jung.Rdata"

vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )

filpath_detr <- c(  paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s0_DETR.nc" ), 
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_DETR.nc" ),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_DETR.nc" ),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_DETR.nc" )
                    )

filpath_nice <- c(  paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s0_ANN.nc" ), 
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_ANN.nc" ),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_ANN.nc" ),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_ANN.nc" )
                    )

modl <- c( "Pmodel_S0", "Pmodel_S1a", "Pmodel_S1b", "Pmodel_S1c")

if (!file.exists(outfiln) || overwrite){

	ampl_1a <- data.frame()
	ampl_1b <- data.frame()
	ampl_1c <- data.frame()

	for ( ires in 1:length(vec_res) ){

	  print(paste("resolution number", ires))
	  
		detr <- list()
		nice <- list()
		for (idx in seq(length(modl))){

		  ## Read detrended GPP for IAV
		  if (file.exists(filpath_detr[idx])){

		    ## read file detrended file
		    if (ires==1){
			    filn <- filpath_detr[idx]
			  } else {
			    filn <- gsub( "DETR", paste0("DETR_REGR", sprintf("%02d", ires)), filpath_detr[idx] )
		    }
		    nc <- nc_open( filn )
		    detr[[ modl[idx] ]] <- ncvar_get( nc, varid="gpp" )
		    nc_close(nc)

		    ## read file non-detrended file
		    if (ires==1){
			    filn <- filpath_nice[idx]
			  } else {
			    filn <- gsub( "ANN", paste0("ANN_REGR", sprintf("%02d", ires)), filpath_nice[idx] )
		    }
		    nc <- nc_open( filn )
		    nice[[ modl[idx] ]] <- ncvar_get( nc, varid="gpp" )
		    nc_close(nc)    

		  }
		  
		}

		## Get variance and relative variance arrays for all simulations
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
		
		ampl_arr_1a <- relvar[[ "Pmodel_S1a" ]] / relvar[[ "Pmodel_S0" ]]
		ampl_arr_1b <- relvar[[ "Pmodel_S1b" ]] / relvar[[ "Pmodel_S0" ]]
		ampl_arr_1c <- relvar[[ "Pmodel_S1c" ]] / relvar[[ "Pmodel_S0" ]]

		tmp_relvar_1a <- c( ampl_arr_1a )
		tmp_relvar_1b <- c( ampl_arr_1b )
		tmp_relvar_1c <- c( ampl_arr_1c )

		idxs <- which( !is.na(tmp_relvar_1a) | !is.na(tmp_relvar_1b) | !is.na(tmp_relvar_1c) )

		tmp_relvar_1a <- tmp_relvar_1a[ idxs ]
		tmp_relvar_1b <- tmp_relvar_1b[ idxs ]
		tmp_relvar_1c <- tmp_relvar_1c[ idxs ]

		addrows <- data.frame( relvar=tmp_relvar_1a, resnr=rep(ires, length(idxs)) )
		ampl_1a <- rbind( ampl_1a, addrows )	

		addrows <- data.frame( relvar=tmp_relvar_1b, resnr=rep(ires, length(idxs)) )
		ampl_1b <- rbind( ampl_1b, addrows )	

		addrows <- data.frame( relvar=tmp_relvar_1c, resnr=rep(ires, length(idxs)) )
		ampl_1c <- rbind( ampl_1c, addrows )	

	}

	save( ampl_1a, ampl_1b, ampl_1c, file=outfiln )

} else {

	load( outfiln )

}

## aggregate to get mean (median and quantiles) amplification
ampl_agg_1a <- ampl_1a %>%  group_by( resnr ) %>%
														summarise(  
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
ampl_agg_1b <- ampl_1b %>%  group_by( resnr ) %>%
														summarise(  
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
ampl_agg_1c <- ampl_1c %>%  group_by( resnr ) %>%
														summarise(  
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

# ## take mean across ensembles
# ampl_agg <- abind( ampl_agg_1a, ampl_agg_1b, ampl_agg_1c, along = 3 ) %>%
# 						apply( c(1,2), FUN = mean ) %>%
#             as_tibble()

## use s1b as standard
ampl_agg <- ampl_agg_1b

save( ampl_agg, file="data/ampl_agg_jung.Rdata" )

## Plot scale dependence of soil moisture effect on GPP interannual variance
pdf("fig/plot_jung.pdf", width = 5, height = 4 )
par(las=1, mar=c(4,4,1,1))
with( ampl_agg, plot( resnr, relvar_median, type="l", col="black", lwd=2, ylim=c(0,10), xaxt = "n", xlab="spatial resolution (degrees)", ylab="ampl. of rel. var. of ann. GPP" ) )
axis( 1, at=seq(length(vec_res)), labels=as.character( vec_res ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q01, relvar_q99), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q05, relvar_q95), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q10, relvar_q90), col=add_alpha("tomato2", 0.25), border = NA ) )
with( ampl_agg, polygon( c(rev(resnr), resnr), c(relvar_q25, relvar_q75), col=add_alpha("tomato2", 0.25), border = NA ) )
abline( h = 1.0, lty=3 )
dev.off()




