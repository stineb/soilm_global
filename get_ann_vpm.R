source("~/.Rprofile")
library(raster)
library(dplyr)
library(abind)

yearstart = 2000
yearend   = 2016

for (year in yearstart:yearend){

  rasta <- raster( paste0( myhome, "data/gpp_vpm/GPP.VPM.", as.character(year), ".v20.HD.tif" ) )
  arr <-  rasta %>% 
            as.array() %>% aperm( perm = c(2,1,3) )
  arr <- arr[,,1]
  
  if (year==2000){
    gpp <- arr
    lon <- seq(-179.75, 179.75, by=0.5)
    lat <- seq(-89.75,89.75, by=0.5)
  } else {
    gpp <- abind( gpp, arr, along=3 )
  }

}

## write annual GPP to file
outfilnam <- paste0( myhome, "/data/gpp_vpm/gpp_vpm_ann.nc")
cdf.write( gpp, "gpp", 
           lon, rev(lat),
           filnam = outfilnam,
           nvars = 1,
           time = yearstart:yearend,
           make.tdim = TRUE,
           units_time = "year",
           long_name_var1 = "Gross primary productivity",
           units_var1 = "gC m-2 year-1",
           glob_hist = "created using soilm_global/get_ann_vpm.R based on original files data/gpp_vpm/GPP.VPM.*.v20.HD.tif, downloaded from https://figshare.com/collections/A_global_moderate_resolution_dataset_of_gross_primary_production_of_vegetation_for_2000-2016/3789814."
           )
