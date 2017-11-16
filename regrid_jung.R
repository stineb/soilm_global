vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )

filpath_detr <- c(  "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_detr_s0_fapar3g_global.nc", 
                    "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_detr_s1_fapar3g_global.nc"
                    )

filpath_nice <- c(  "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_nice_s0_fapar3g_global.nc", 
                    "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_nice_s1_fapar3g_global.nc"
                    )

## regrid using 'cdo remapbil'
for ( ires in 2:(length(vec_res)-1) ){

	for (isim in 1:2){

		## regrid detrended file
		gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
		infile <- filpath_detr[isim]
		outfile <- gsub( "detr", paste0("detr_regr", sprintf("%02d", ires)), infile )
		cmd <- paste0( "cdo remapbil,", gridfile, " ", infile, " ", outfile )
		print( cmd )
		system( cmd )

		## regrid nice file
		gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
		infile <- filpath_nice[isim]
		outfile <- gsub( "nice", paste0("nice_regr", sprintf("%02d", ires)), infile )
		cmd <- paste0( "cdo remapbil,", gridfile, " ", infile, " ", outfile )
		print( cmd )
		system( cmd )	

	}

}

## the last one (global) has to be treated a bit differently using 'cdo fldmean'
ires <- length(vec_res)
for (isim in 1:2){

	## regrid detrended file
	gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
	infile <- filpath_detr[isim]
	outfile <- gsub( "detr", paste0("detr_regr", sprintf("%02d", ires)), infile )
	cmd <- paste0( "cdo fldmean ", infile, " ", outfile )
	print( cmd )
	system( cmd )

	## regrid nice file
	gridfile <- paste0("grids/grid_", sprintf("%02d", ires),".txt")
	infile <- filpath_nice[isim]
	outfile <- gsub( "nice", paste0("nice_regr", sprintf("%02d", ires)), infile )
	cmd <- paste0( "cdo fldmean ", infile, " ", outfile )
	print( cmd )
	system( cmd )	

}
