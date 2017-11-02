complment_all <- function( linearfit, nlsfit ){

	source("compl_df_flue_est.R")

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'
	nice_agg <- compl_df_flue_est( nice_agg, linearfit, nlsfit )
	save( nice_agg, file="data/nice_nn_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, nlsfit )
	save( modis_agg, file="data/nice_nn_modis_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_modis_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, nlsfit )
	save( modis_agg, file="data/nice_nn_modis_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'
	nice_agg <- compl_df_flue_est( nice_agg, linearfit, nlsfit )
	save( nice_agg, file="data/nice_all_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_all_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, nlsfit )
	save( modis_agg, file="data/nice_all_modis_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_all_modis_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, nlsfit )
	save( modis_agg, file="data/nice_all_modis_agg_lue_obs_evi_L2.Rdata" )


	## complement nice files
	files <- list.files("data", pattern = "nice_nn_*") 
	files <- files[-grep("L2", files, fixed=T)]
	files <- files[-grep("agg_lue_obs_evi", files, fixed=T)]

	files2 <- list.files("data", pattern = "nice_all_*") 
	files2 <- files2[-grep("L2", files2, fixed=T)]
	files2 <- files2[-grep("agg_lue_obs_evi", files2, fixed=T)]

	files <- c( files, files2 )

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice <- compl_df_flue_est( nice, linearfit, nlsfit )
		save( nice, file=paste0( "data/", ifil ) )

	}


	## complement modis files
	files <- list.files("data", pattern = "modis_*") 
	files <- files[-grep("df_dday", files, fixed=T)]
	files <- files[-grep("nice", files, fixed=T)]

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice_to_modis <- compl_df_flue_est( nice_to_modis, linearfit, nlsfit )
		save( nice_to_modis, file=paste0( "data/", ifil ) )

	}


	## complement mte files
	files <- list.files("data", pattern = "mte_*") 
	files <- files[-grep("df_dday", files, fixed=T)]
	files <- files[-grep("nice", files, fixed=T)]

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice_to_mte <- compl_df_flue_est( nice_to_mte, linearfit, nlsfit )
		save( nice_to_mte, file=paste0( "data/", ifil ) )

	}

}