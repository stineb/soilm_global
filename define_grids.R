
vec_res <- c( 0.5, 1.0, 1.5, 2.5, 3, 4, 4.5, 5, 6, 7.5, 9, 10, 12, 15, 18, 20, 22.5, 30, 36, 45, 60, 90, 180, 360 )

for (idx in seq(length(vec_res))){

	dlat <- vec_res[idx]
	dlon <- dlat
	nlon <- 360 / dlon
	nlat <- 180 / dlat

	con <- file( paste0("./grids/grid_", sprintf("%02d", idx),".txt"), "w" )
	cat( "gridtype = ", "lonlat",    "\n", file=con, sep=" ")
	cat( "xsize    = ", nlon,        "\n", file=con, sep=" ")
	cat( "ysize    = ", nlat,        "\n", file=con, sep=" ")
	cat( "xfirst   = ", -180+dlon/2, "\n", file=con, sep=" ")
	cat( "xinc     = ", dlon,        "\n", file=con, sep=" ")
	cat( "yfirst   = ", -90+dlat/2,  "\n", file=con, sep=" ")
	cat( "yinc     = ", dlat,        "\n", file=con, sep=" ")
	close(con)

}
