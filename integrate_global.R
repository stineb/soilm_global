integrate_global <- function( arr ){

  if ( dim(arr)[1]==720 && dim(arr)[2]==360 ){

    ## half degree resolution
    nc <- nc_open( "data/area_halfdeg.nc" )
    arr_area <- ncvar_get( nc, varid="area" )
    nc_close( nc )  

  } else if ( dim(arr)[1]==360 && dim(arr)[2]==180 ){

    ## one degree resolution
    ## half degree resolution
    nc <- nc_open( "data/area_1x1deg.nc" )
    arr_area <- ncvar_get( nc, varid="area" )
    nc_close( nc )  

  }    

  if (length(dim(arr))==2){
    
    arr_abs <- arr * arr_area
    garr <- sum( arr_abs, na.rm=TRUE )

  } else if (length(dim(arr))==3){

    ## get global total unlimited GPP over time
    arr_abs <- sweep( arr, 1, arr_area, "*", check.margin=FALSE )
    garr <- apply( arr_abs, c(3), FUN=sum, na.rm=TRUE )

  } else {
    print("cannot deal with this number of dimensions")
    garr <- NA
  }

  return( garr )

}