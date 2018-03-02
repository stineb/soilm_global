library(ncdf4)
library(abind)

source("get_ahlstroem_f.R")
source("../utilities/plot_map.R")

filpath_detr <- c(  paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s0_DETR.nc"), 
                    paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1a_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1b_DETR.nc"),
                    paste0( myhome, "/data/pmodel_fortran_output/gpp_pmodel_s1c_DETR.nc")
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
ahlstroem_fa <- get_ahlstroem_f( detr$diffa, isabs=FALSE )
ahlstroem_fb <- get_ahlstroem_f( detr$diffb, isabs=FALSE )
ahlstroem_fc <- get_ahlstroem_f( detr$diffc, isabs=FALSE )

ahlstroem_f <-  abind( ahlstroem_fa, ahlstroem_fb, ahlstroem_fc, along = 3 ) %>%
				        apply( c(1,2), FUN = mean ) 

par( mgp=c(3,1,0) )

plot_map( ahlstroem_fa*1e4, lev=seq(-4.5,4.5,0.5), positive=FALSE, maxval=4 ) # , file="fig/map_ahlstroem_gpploss.pdf"

