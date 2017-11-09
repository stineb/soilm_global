complment_all <- function( linearfit ){

	source("compl_df_flue_est.R")

	x0_fix = 0.9

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'
	nice_agg <- compl_df_flue_est( nice_agg, linearfit, x0_fix=x0_fix )
	save( nice_agg, file="data/nice_nn_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, x0_fix=x0_fix )
	save( modis_agg, file="data/nice_nn_modis_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_modis_agg_lue_obs_evi.Rdata'
	load( "data/nice_nn_mte_agg_lue_obs_evi.Rdata" )  # loads 'mte_agg'
	mte_agg <- compl_df_flue_est( mte_agg, linearfit, x0_fix=x0_fix )
	save( mte_agg, file="data/nice_nn_mte_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_nn_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_agg_lue_obs_evi.Rdata" )  # loads 'nice_agg'
	nice_agg <- compl_df_flue_est( nice_agg, linearfit, x0_fix=x0_fix )
	save( nice_agg, file="data/nice_all_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_all_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_modis_agg_lue_obs_evi.Rdata" )  # loads 'modis_agg'
	modis_agg <- compl_df_flue_est( modis_agg, linearfit, x0_fix=x0_fix )
	save( modis_agg, file="data/nice_all_modis_agg_lue_obs_evi_L2.Rdata" )

	## complementing 'data/nice_all_modis_agg_lue_obs_evi.Rdata'
	load( "data/nice_all_mte_agg_lue_obs_evi.Rdata" )  # loads 'mte_agg'
	mte_agg <- compl_df_flue_est( mte_agg, linearfit, x0_fix=x0_fix )
	save( mte_agg, file="data/nice_all_mte_agg_lue_obs_evi_L2.Rdata" )


	# ## complementing 'data/nice_all_modis_agg_lue_obs_evi.Rdata'
	# load( "data/data_aligned_agg.Rdata" )  # loads 'mte_agg'
	# df_dday_agg <- compl_df_flue_est( df_dday_agg, linearfit, x0_fix=x0_fix )
	# df_dday_modis_agg <- compl_df_flue_est( df_dday_modis_agg, linearfit, x0_fix=x0_fix )
	# df_dday_mte_agg <- compl_df_flue_est( df_dday_mte_agg, linearfit, x0_fix=x0_fix )
	# save( df_dday_agg, df_dday_modis_agg, df_dday_mte_agg, df_dday_aggbydday_agg, file="data/data_aligned_agg_L2.Rdata" )


	## complement nice files
	files <- list.files("data", pattern = "nice_nn_*") 
	files <- files[-grep("L2", files, fixed=TRUE )]
	files <- files[-grep("agg_lue_obs_evi", files, fixed=TRUE )]

	files2 <- list.files("data", pattern = "nice_all_*") 
	files2 <- files2[-grep("L2", files2, fixed=TRUE )]
	files2 <- files2[-grep("agg_lue_obs_evi", files2, fixed=TRUE )]

	files <- c( files, files2 )

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice <- compl_df_flue_est( nice, linearfit, x0_fix=x0_fix )
		save( nice, file=paste0( "data/", ifil ) )

	}

	## complement modis files
	files <- list.files("data", pattern = "modis_*") 
	files <- files[-grep("df_dday", files, fixed=TRUE )]
	files <- files[-grep("nice", files, fixed=TRUE )]

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice_to_modis <- compl_df_flue_est( nice_to_modis, linearfit, x0_fix=x0_fix )
		save( nice_to_modis, file=paste0( "data/", ifil ) )

	}


	## complement mte files
	files <- list.files("data", pattern = "mte_*") 
	files <- files[-grep("df_dday", files, fixed=TRUE )]
	files <- files[-grep("nice", files, fixed=TRUE )]

	for (ifil in files){

		load( paste0( "data/", ifil ) )
		nice_to_mte <- compl_df_flue_est( nice_to_mte, linearfit, x0_fix=x0_fix )
		save( nice_to_mte, file=paste0( "data/", ifil ) )

	}

}