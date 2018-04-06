vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )

filpath_DETR <- c(  "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s0_DETR.nc", 
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_DETR.nc",
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_DETR.nc",
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_DETR.nc"
                    )

filpath_ANN <- c(  "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s0_ANN.nc", 
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_ANN.nc",
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_ANN.nc",
                    "/Users/benjaminstocker/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_ANN.nc"
                    )

## regrid using 'cdo remapbil'
for ( ires in 2:(length(vec_res)-1) ){

	for (isim in 1:length(filpath_DETR)){

		## regrid detrended file
		gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
		infile <- filpath_DETR[isim]
		outfile <- gsub( "DETR", paste0("DETR_REGR", sprintf("%02d", ires)), infile )
		cmd <- paste0( "cdo remapbil,", gridfile, " ", infile, " ", outfile )
		print( cmd )
		system( cmd )

		## regrid ANN file
		gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
		infile <- filpath_ANN[isim]
		outfile <- gsub( "ANN", paste0("ANN_REGR", sprintf("%02d", ires)), infile )
		cmd <- paste0( "cdo remapbil,", gridfile, " ", infile, " ", outfile )
		print( cmd )
		system( cmd )	

	}

}

## the last one (global) has to be treated a bit differently using 'cdo fldmean'
ires <- length(vec_res)
for (isim in 1:length(filpath_DETR)){

	## regrid detrended file
	gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
	infile <- filpath_DETR[isim]
	outfile <- gsub( "DETR", paste0("DETR_REGR", sprintf("%02d", ires)), infile )
	cmd <- paste0( "cdo fldmean ", infile, " ", outfile )
	print( cmd )
	system( cmd )

	## regrid ANN file
	gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
	infile <- filpath_ANN[isim]
	outfile <- gsub( "ANN", paste0("ANN_REGR", sprintf("%02d", ires)), infile )
	cmd <- paste0( "cdo fldmean ", infile, " ", outfile )
	print( cmd )
	system( cmd )	

}
