get_ahlstroem_f <- function( arr, isabs=FALSE ){
  ##---------------------------------------------------
  ## requires as input a 3D array with lon x lat x time
  ## and values being some interannual variability 
  ## index quantifies the degree to which each gridcell contributes to the global signal
  ##---------------------------------------------------

  ## use gridcell totals not per unit area
  if (isabs==FALSE){
    source( "integrate_gridcell.R" )
    arr_abs <- integrate_gridcell( arr, global=FALSE, overwrite=FALSE )
    glob <- apply( arr_abs, c(3), FUN=sum, na.rm=TRUE )
  } else {
    arr_abs <- arr
    glob <- apply( arr_abs, c(3), FUN=sum, na.rm=TRUE )    
  }

  ahlstroem_f <- arr[,,1]
  ahlstroem_f[] <- NA
  for (ilon in seq(dim(arr)[1])){
    for (ilat in seq(dim(arr)[2])){
      if (!is.na(arr[ilon,ilat,1])){
        ahlstroem_f[ilon,ilat] <- sum( arr_abs[ilon,ilat,] * abs( glob ) / glob ) / sum( abs( glob ) )
      }
    }
  }
  return( ahlstroem_f )
}