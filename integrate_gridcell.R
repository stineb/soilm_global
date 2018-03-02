integrate_gridcell <- function( arr, global=TRUE, overwrite=FALSE ){

  ## get area array
  if ( dim(arr)[1]==720 && dim(arr)[2]==360 ){

    ## half degree resolution
    print("found halfdegree resolution.")
    areafil <- paste0( myhome, "data/landmasks/area_halfdeg.nc")

    if (file.exists(areafil)&&!overwrite){
   
      print("reading from file.")     
      nc <- nc_open( areafil )
      arr_area <- ncvar_get( nc, varid="area" )
      nc_close( nc )        

    } else {
      
      dx <- 0.5
      dy <- 0.5
      lon <- seq(-179.75, 179.75, by=dx )
      lat <- seq(-89.75, 89.75, by=dy )
      if (length(dim(arr))==3){
        arr_area <- arr[,,1]
        arr_tmp  <- arr[,,1]
      } else if (length(dim(arr))==2) {
        arr_area <- arr
        arr_tmp  <- arr
      } else {
        print("cannot deal with this number of dimensions")
      }
      arr_area[] <- NA
      for (ilon in seq(dim(arr)[1])){
        for (ilat in seq(dim(arr)[2])){
          if (!is.na(arr_tmp[ilon,ilat])){
            arr_area[ilon,ilat] <- area( lat[ilat], dx=dx, dy=dy )
          }
        }
      }
      cdf.write( arr_area, "area", 
                 lon, lat,
                 filnam = areafil,
                 nvars = 1,
                 make.tdim = FALSE,
                 long_name_var1 = "gridcell area",
                 units_var1 = "m2",
                 glob_hist = "created by soilm_global/integrate_global.R",
                 glob_title = "gridcell area"
      )
    }

  } else if ( dim(arr)[1]==360 && dim(arr)[2]==180 ){

    ## one degree resolution
    print("found one degree resolution.")
    areafil <- paste0( myhome, "/data/landmasks/area_1x1deg.nc")

    if (file.exists(areafil)&&!overwrite){

      nc <- nc_open( areafil )
      arr_area <- ncvar_get( nc, varid="area" )
      nc_close( nc )        

    } else {

      dx <- 1.0
      dy <- 1.0
      lon <- seq(-175.5, 175.5, dx )
      lat <- seq(-89.5, 89.5, dy )
      arr_area <- arr[,,1]
      arr_area[] <- NA
      for (ilon in seq(dim(arr)[1])){
        for (ilat in seq(dim(arr)[2])){
          if (!is.na(arr[ilon,ilat,])){
            arr_area[ilon,ilat] <- area( lat[ilat], dx=dx, dy=dy )
          }
        }
      }
      cdf.write( arr_area, "area", 
                 lon, lat,
                 filnam = areafil,
                 nvars = 1,
                 make.tdim = FALSE,
                 long_name_var1 = "gridcell area",
                 units_var1 = "m2",
                 glob_hist = "created by soilm_global/integrate_global.R",
                 glob_title = "gridcell area"
      )
    }

  }

  if (!global){

    ## actually integrate
    if (length(dim(arr))==2){
      
      ## 2D array    
      out <- arr * arr_area

    } else if (length(dim(arr))==3){
      
      ## 3D arry
      ## get global total unlimited GPP over time
      out <- sweep( arr, c(1,2), arr_area, "*", check.margin=FALSE )

    } else {

      print("cannot deal with this number of dimensions")
      out <- NA

    }
    
  } else {

    ## actually integrate
    if (length(dim(arr))==2){
      
      ## 2D array
      arr_abs <- arr * arr_area
      out <- sum( arr_abs, na.rm=TRUE )

    } else if (length(dim(arr))==3){

      ## 3D arry
      ## get global total unlimited GPP over time
      arr_abs <- sweep( arr, c(1,2), arr_area, "*", check.margin=FALSE )
      out <- apply( arr_abs, c(3), FUN=sum, na.rm=TRUE )
  
      # ## This is a correct application of sweep, as demonstrated below...
      # arr_abs_test <- arr * NA
      # for (itim in seq(dim(arr)[3])){
      #   arr_abs_test[,,itim] <- arr[,,itim] * arr_area[,]
      # }
      # print(all.equal( arr_abs, arr_abs_test ))
      
      
    } else {

      print("cannot deal with this number of dimensions")
      out <- NA

    }

  }

  return( out )

}