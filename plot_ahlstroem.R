library(ncdf4)
library(abind)

source("get_ahlstroem_f.R")
source("../utilities/plot_map.R")


get_stocker_f <- function( eff, anom, isabs=FALSE ){
  ##---------------------------------------------------
  ## requires as input a 3D effay with lon x lat x time
  ## and values being annual detrended anomalies 
  ## index quantifies the degree to which each gridcell contributes to the global signal
  ##---------------------------------------------------

  ## use gridcell total not per unit area
  if (isabs==FALSE){
    source( "integrate_gridcell.R" )
    
    eff_abs <- integrate_gridcell( eff, global=FALSE, overwrite=TRUE )
    anom_abs <- integrate_gridcell( anom, global=FALSE, overwrite=TRUE )

    eff_glob <- apply( eff_abs, c(3), FUN=sum, na.rm=TRUE )
    anom_glob <- apply( anom_abs, c(3), FUN=sum, na.rm=TRUE )
  } else {
    
    eff_abs <- eff
    eff_glob <- apply( eff_abs, c(3), FUN=sum, na.rm=TRUE )    
    anom_glob <- apply( anom_abs, c(3), FUN=sum, na.rm=TRUE )    
  }

  stocker_f <- eff[,,1]
  stocker_f[] <- NA
  for (ilon in seq(dim(eff)[1])){
    for (ilat in seq(dim(eff)[2])){
      if (!is.na(eff[ilon,ilat,1])){
        # ## version 1
        # stocker_f[ilon,ilat] <- sum( eff_abs[ilon,ilat,] * abs( anom_glob ) / eff_glob ) / sum( abs( anom_glob ) )

        ## version 2
        stocker_f[ilon,ilat] <- sum( eff_abs[ilon,ilat,] * abs( anom_glob ) / anom_glob ) / sum( abs( anom_glob ) )
      }
    }
  }
  return( stocker_f )
}


filpath_detr <- c(  paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s0_DETR.nc"), 
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_DETR.nc")
                    )

modl <- c( "Pmodel_S0", "Pmodel_S1a", "Pmodel_S1b", "Pmodel_S1c")

detr <- list()
for (idx in seq(length(modl))){

  ## Read detrended GPP for IAV
  if (file.exists(filpath_detr[idx])){

    ## read file
    nc <- nc_open( filpath_detr[idx] )
    detr[[ modl[idx] ]] <- try( ncvar_get( nc, varid="gpp" ) )
    nc_close(nc)

  }
  
}

## get difference
detr[[ "diffa" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1a" ]]
detr[[ "diffb" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1b" ]]
detr[[ "diffc" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1c" ]]

## get ahlstroem-f
# ahlstroem_fa <- get_ahlstroem_f( detr$diffa, isabs=FALSE )
ahlstroem_fb <- get_ahlstroem_f( detr$diffb, isabs=FALSE )
# ahlstroem_fc <- get_ahlstroem_f( detr$diffc, isabs=FALSE )

ahlstroem_f <-  abind( ahlstroem_fa, ahlstroem_fb, ahlstroem_fc, along = 3 ) %>%
				        apply( c(1,2), FUN = mean ) 

par( mgp=c(3,1,0) )

plot_map( ahlstroem_fb*1e4, lev=seq(-2.5,2.5,0.5), positive=FALSE, maxval=4, file="fig/map_ahlstroem_gpploss.pdf" ) #

# ## version 1
# stocker_fb <- get_stocker_f( detr$diffb, detr$Pmodel_S0, isabs=FALSE )
# plot_map( stocker_fb*1e4, lev=seq(-5,5,0.5), positive=FALSE, maxval=30, minval=-30, file="fig/map_stocker_gpploss.pdf" ) #

## version 2
stocker_fb <- get_stocker_f( detr$diffb, detr$Pmodel_S1b, isabs=FALSE )
plot_map( stocker_fb*1e4, lev=seq(-0.5,0.5,0.1), positive=FALSE, maxval=30, minval=-30, file="fig/map_stocker_gpploss.pdf" ) #



