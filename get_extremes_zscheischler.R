# make anomalies
require(raster)
require(neuroim)
require(dplyr)
require(ncdf4)
require(lubridate)

myboxplot <- function( ... ){
  boxplot( ..., staplewex=0, whisklty=1, outline=FALSE )
}


##-----------------------
overwrite <- TRUE
##-----------------------

outfile <- "data/extremes.Rdata"
nconts <- 5

##------------------------------------------------------------
## Identify events
##------------------------------------------------------------
if (!file.exists(outfile)||overwrite){

  ## Define continents from SREX regions
  # SREX <- raster("/net/exo/landclim/data/dataset/SREX-Region-Masks/20120709/0.5deg_lat-lon_time-invariant/processed/netcdf/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
  SREX <- raster("~/data/landmasks/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
  subsets <- list()

  # subsets[[1]] <- 1:6    # north plus central america
  # subsets[[2]] <- 7:10   # south america
  # subsets[[3]] <- 11:13  # europe
  # subsets[[4]] <- 14:17  # africa
  # subsets[[5]] <- 18:23  # asia excluding southeast asia (24)
  # subsets[[6]] <- 24:26  # australia

  subsets[[1]] <- 1:6    # north plus central america
  subsets[[2]] <- 7:10   # south america
  subsets[[3]] <- c(11:13, 18:23)  # europe + asia
  subsets[[4]] <- 14:17  # africa
  subsets[[5]] <- 24:26  # australia

  nconts <- length(subsets)

  CC_list <- list()

  idx_largest <- list()

  fit_s0  <- list()
  fit_s1a <- list()
  fit_s1b <- list()
  fit_s1c <- list()

  IMPACT_s0  <- list()
  IMPACT_s1a <- list()
  IMPACT_s1b <- list()
  IMPACT_s1c <- list()

  ## load monthly anomalies created beforehand with 'proc_nc_fields.sh'
  print("loading anomaly files...")
  ANOM_s0_all  <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s0_MON_ANOM.nc")
  ANOM_s1a_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_MON_ANOM.nc")
  ANOM_s1b_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_MON_ANOM.nc")
  ANOM_s1c_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_MON_ANOM.nc")

  ## reorder dimensions and multiply with gridcell surface area
  print("reorder dimensions...")
  ANOM_s0.a  <- aperm(as.array(ANOM_s0_all  * area(ANOM_s0_all)),c(2,1,3))[,360:1,]
  ANOM_s1a.a <- aperm(as.array(ANOM_s1a_all * area(ANOM_s1a_all)),c(2,1,3))[,360:1,]
  ANOM_s1b.a <- aperm(as.array(ANOM_s1b_all * area(ANOM_s1b_all)),c(2,1,3))[,360:1,]
  ANOM_s1c.a <- aperm(as.array(ANOM_s1c_all * area(ANOM_s1c_all)),c(2,1,3))[,360:1,]

  ## Loop over continents
  latx <- 360
  for (j in 1:nconts) {
    
    print(paste("continent",as.character(j),"/ ", as.character(nconts) ))
    
    # define mask for only 1 continent
    Na <- SREX %in% subsets[[j]]

    ## mask out values
    ANOM_s1a <- mask( ANOM_s1a_all, Na, maskvalue=0 )
    ANOM_s1b <- mask( ANOM_s1b_all, Na, maskvalue=0 )
    ANOM_s1c <- mask( ANOM_s1c_all, Na, maskvalue=0 )
    
    ## get quantiles based on simulation s1c
    q_all <- quantile( as.vector(ANOM_s1b), c(0.01,0.02,0.03,0.05,0.1,0.9,0.95), na.rm = TRUE )
    
    # same quantiles (computed over the combined dataset)
    EXT05_s1b <- ANOM_s1b < q_all[2]
    
    EXT05_s1b <- aperm( as.array(EXT05_s1b), c(2,1,3) )[,latx:1,] == TRUE
    print("       identifying events...")
    CC_s1 <- connComp3D( EXT05_s1b )
    
    ## keep only events gt 200
    print("       calculating impacts...")
    CC_s1$index[CC_s1$size < 50] <- 0
    
    # make list with indices pointing to equal index
    # CC_s1$list are the indeces pointing to respective events, it's a list with length equal number of events
    CC_s1$list <- sapply( 1:max(CC_s1$index), function(x) which(CC_s1$index == x) )
    
    # find largest extremes defined by impact
    IMPACT_s0[[j]]  <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s0.a[x])) )
    IMPACT_s1a[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1a.a[x])) )
    IMPACT_s1b[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1b.a[x])) )
    IMPACT_s1c[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1c.a[x])) )

    # attach to list
    CC_list[[j]] <- CC_s1

  }

  print("saving...")
  save( CC_list, IMPACT_s0, IMPACT_s1a, IMPACT_s1b, IMPACT_s1c, file=outfile )

} else {

  load(outfile)

}

##------------------------------------------------------------
## Get size of events (impact), dataframe 
##------------------------------------------------------------
## Size difference of the 10 largest extremes
## sort by descending s1c impacts
print("constructing impacts dataframe...")
list_impacts <- list()
for (icont in 1:nconts){
  df_impacts <- tibble( s0 = IMPACT_s0[[icont]], s1a = IMPACT_s1a[[icont]], s1b = IMPACT_s1b[[icont]], s1c = IMPACT_s1c[[icont]] ) %>%
    filter( s0<0 & s1a<0 & s1b<0 & s1c<0 ) %>%
    mutate( rank     = rank( s1b ),
            ampl_s1a = s1a / s0, 
            ampl_s1b = s1b / s0, 
            ampl_s1c = s1c / s0,
            icont    = icont,
            event_no = 1:nrow(.)
          ) %>%
    mutate( ampl_mean = rowMeans( dplyr::select(., starts_with("ampl")) ) )
  list_impacts[[icont]] <- df_impacts
}

## combine list into single dataframe (tibble) and get rank by size across all continents (impact in s1b)
df_impacts <- bind_rows( list_impacts ) %>% 
              mutate( rank_global = rank( s1b ) ) %>%
              arrange( s1b )
              

##------------------------------------------------------------
## Save all data needed for plotting later (plot_fig_5.R)
##------------------------------------------------------------
save( df_impacts, list_impacts, ANOM_s0.a, ANOM_s1b.a, CC_list, file = "data/impacts_extremes.Rdata")


##------------------------------------------------------------
## Locate events in time and space (only 141 biggest)
##------------------------------------------------------------
## integrate over ten largest events globally
## get dates from NetCDF file
print("plot events...")
nc <- nc_open( "~/data/pmodel_fortran_output/v2/gpp_pmodel_s0_MON_ANOM.nc" )
nc_close(nc)
date <- ymd( "2001-01-01" ) + days( floor(nc$dim$time$vals) )
nevents <- 200

overwrite_located <- TRUE
filn_located <- "data/extremes_located.Rdata"

if (file.exists(filn_located)&&!overwrite_located){
  
  load( filn_located )

} else {

  list_event_time   <- list()
  list_event_lonlat <- list()
  
  df_event_time <- tibble()
  arr_event_lonlat  <- array( NA, dim = c(720,360,nevents) )
  addmap <- FALSE

  for (irank in 1:nevents){
    
    print( paste( "    event ranked", as.character(irank) ) )
    tmp <- df_impacts %>% filter( rank_global==irank )
    myicont <- tmp$icont
    myevent_no <- tmp$event_no
    rank_global <- tmp$rank_global
    
    ## create array that contains anomaly for voxels in this event and NA otherwise
    ## ... based on s0
    LargestE_s0 <- NA * ANOM_s0.a
    LargestE_s0[CC_list[[myicont]]$list[[myevent_no]]] <- ANOM_s0.a[CC_list[[myicont]]$list[[myevent_no]]]
    
    ## ... based on s1b
    LargestE_s1b <- NA * ANOM_s1b.a
    LargestE_s1b[CC_list[[myicont]]$list[[myevent_no]]] <- ANOM_s1b.a[CC_list[[myicont]]$list[[myevent_no]]]
    
    ## integrate over lon-lat, projecting onto time
    event_time_s0    <- apply( LargestE_s0 , 3, sum, na.rm = TRUE ) * 1e-15
    event_time_s1b   <- apply( LargestE_s1b, 3, sum, na.rm = TRUE ) * 1e-15

    addrows <- tibble( date=date, year=year(date), month=month(date), impact_s0=event_time_s0, impact_s1b=event_time_s1b, icont=myicont, event_no=myevent_no, rank_global=rank_global ) %>%
               filter( abs(impact_s0)>0 | abs(impact_s1b)>0 )
    df_event_time <- bind_rows( df_event_time, addrows )

    ## integrate over time, projecting into lon-lat, based on s1b
    arr_event_lonlat[,,irank] <- apply( LargestE_s1b, 1:2, sum, na.rm = TRUE )
    arr_event_lonlat[ arr_event_lonlat==0 ] <- NA

    ## plot map of impact, integrated over time
    yearchar <- df_event_time %>% dplyr::filter( event_no==myevent_no & icont==myicont ) %>% mutate( year=year(date)) %>% dplyr::select(year)  %>% table() %>% sort() %>% names()
    if (irank <= nevents_plot) plot_map( -arr_event_lonlat[,,irank], col = cols[irank], add=addmap, file=paste0( "fig/map_event_", sprintf( "%02d", irank ), ".pdf" ), text=yearchar[1] )
    
  }
  
  save( arr_event_lonlat, df_event_time, file=filn_located )
  
}
