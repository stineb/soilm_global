library(ncdf4)

source("get_ahlstroem_f.R")
source("../utilities/plot_map.R")

filpath_detr <- c(  "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_detr_s0_fapar3g_global.nc", 
                    "/Users/benjaminstocker/data/pmodel_fortran_output/pmodel_gpp_detr_s1_fapar3g_global.nc"
                    )

modl <- c( "Pmodel_S0", "Pmodel_S1")

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
detr[[ "diff" ]] <- detr[[ "Pmodel_S0" ]] - detr[[ "Pmodel_S1" ]]

## get ahlstroem-f
ahlstroem_f <- get_ahlstroem_f( detr$diff, isabs=FALSE )

plot_map( ahlstroem_f*1e4, lev=seq(-2,2,0.25), positive=FALSE, file="fig/map_ahlstroem_gpploss.pdf" )
