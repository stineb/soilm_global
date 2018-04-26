# # make anomalies
# require(raster)
# require(neuroim)
# require(igraph)
# require(fields) ## for image.plot() function
# require(dplyr)
# require(ncdf4)
# require(lubridate)
# require(fields)
# require(sp)
# require(maptools)
# 
# myboxplot <- function( ... ){
#   boxplot( ..., staplewex=0, whisklty=1, outpch=16, outcex=0.5 )
# }
# 
# plot_map <- function( out, col, add, file=NA, text="" ){
#   ylim <- c(-60,85)
#   lon <- seq(-179.75, 179.75, 0.5)
#   lat  <- seq(-89.75, 89.75, 0.5)
#   lat.labels <- seq(-90, 90, 30)
#   lat.short  <- seq(-90, 90, 10)
#   lon.labels <- seq(-180, 180, 60)
#   lon.short  <- seq(-180, 180, 10)
#   a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
#   b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )
#   
#   if (!is.na(file)) pdf( file, width=8, height=6 )
#     par( mar=c(3,3,3,1),xaxs="i", yaxs="i",las=1)
#     image(
#       lon, lat, 
#       out,
#       ylim=c(-60,85), 
#       yaxt="n", xaxt="n",
#       col=colorRampPalette( c( add_alpha(col, 0.3), col ) )(10),
#       xlab="", ylab="",
#       add=add
#     )
#     map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
#     axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
#     axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )
#     axis( 4, at=lat.labels, lab=F, lwd=1.5 )
#     axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )
#     axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
#     axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )
#     axis( 3, at=lon.labels, lab=F, lwd=1.5 )
#     axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )
#     text( -160, -50, text, col=col )
#   if (!is.na(file)) dev.off()  
# 
# }
# 
# ##-----------------------
# overwrite <- FALSE
# ##-----------------------
# 
# outfile <- " "
# 
# ##------------------------------------------------------------
# ## Identify events
# ##------------------------------------------------------------
# if (!file.exists(outfile)||overwrite){
# 
#   ## Define continents from SREX regions
#   # SREX <- raster("/net/exo/landclim/data/dataset/SREX-Region-Masks/20120709/0.5deg_lat-lon_time-invariant/processed/netcdf/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
#   SREX <- raster("~/data/landmasks/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
#   subsets <- list()
# 
#   # subsets[[1]] <- 1:6    # north plus central america
#   # subsets[[2]] <- 7:10   # south america
#   # subsets[[3]] <- 11:13  # europe
#   # subsets[[4]] <- 14:17  # africa
#   # subsets[[5]] <- 18:23  # asia excluding southeast asia (24)
#   # subsets[[6]] <- 24:26  # australia
# 
#   subsets[[1]] <- 1:6    # north plus central america
#   subsets[[2]] <- 7:10   # south america
#   subsets[[3]] <- c(11:13, 18:23)  # europe + asia
#   subsets[[4]] <- 14:17  # africa
#   subsets[[5]] <- 24:26  # australia
# 
#   nconts <- length(subsets)
# 
#   CC_list <- list()
# 
#   idx_largest <- list()
# 
#   fit_s0  <- list()
#   fit_s1a <- list()
#   fit_s1b <- list()
#   fit_s1c <- list()
# 
#   IMPACT_s0  <- list()
#   IMPACT_s1a <- list()
#   IMPACT_s1b <- list()
#   IMPACT_s1c <- list()
# 
#   ## load monthly anomalies created beforehand with 'proc_nc_fields.sh'
#   print("loading anomaly files...")
#   ANOM_s0_all  <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s0_MON_ANOM.nc")
#   ANOM_s1a_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1a_MON_ANOM.nc")
#   ANOM_s1b_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1b_MON_ANOM.nc")
#   ANOM_s1c_all <- brick("~/data/pmodel_fortran_output/v2/gpp_pmodel_s1c_MON_ANOM.nc")
# 
#   ## reorder dimensions and multiply with gridcell surface area
#   print("reorder dimensions...")
#   ANOM_s0.a  <- aperm(as.array(ANOM_s0_all  * area(ANOM_s0_all)),c(2,1,3))[,360:1,]
#   ANOM_s1a.a <- aperm(as.array(ANOM_s1a_all * area(ANOM_s1a_all)),c(2,1,3))[,360:1,]
#   ANOM_s1b.a <- aperm(as.array(ANOM_s1b_all * area(ANOM_s1b_all)),c(2,1,3))[,360:1,]
#   ANOM_s1c.a <- aperm(as.array(ANOM_s1c_all * area(ANOM_s1c_all)),c(2,1,3))[,360:1,]
# 
#   ## Loop over continents
#   latx <- 360
#   for (j in 1:nconts) {
#     
#     print(paste("continent",as.character(j),"/ 6"))
#     
#     # define mask for only 1 continent
#     Na <- SREX %in% subsets[[j]]
# 
#     ## mask out values
#     ANOM_s1a <- mask( ANOM_s1a_all, Na, maskvalue=0 )
#     ANOM_s1b <- mask( ANOM_s1b_all, Na, maskvalue=0 )
#     ANOM_s1c <- mask( ANOM_s1c_all, Na, maskvalue=0 )
#     
#     ## get quantiles based on simulation s1c
#     q_all <- quantile( as.vector(ANOM_s1b), c(0.01,0.03,0.05,0.1,0.9,0.95), na.rm = TRUE )
#     
#     # same quantiles (computed over the combined dataset)
#     EXT05_s1b <- ANOM_s1b < q_all[1]
#     
#     EXT05_s1b <- aperm( as.array(EXT05_s1b), c(2,1,3) )[,latx:1,] == TRUE
#     print("       identifying events...")
#     CC_s1 <- connComp3D( EXT05_s1b )
#     
#     ## keep only events gt 200
#     print("       calculating impacts...")
#     CC_s1$index[CC_s1$size < 200] <- 0
#     
#     # make list with indices pointing to equal index
#     # CC_s1$list are the indeces pointing to respective events, it's a list with length equal number of events
#     CC_s1$list <- sapply( 1:max(CC_s1$index), function(x) which(CC_s1$index == x) )
#     
#     # find largest extremes defined by impact
#     IMPACT_s0[[j]]  <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s0.a[x])) )
#     IMPACT_s1a[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1a.a[x])) )
#     IMPACT_s1b[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1b.a[x])) )
#     IMPACT_s1c[[j]] <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1c.a[x])) )
# 
#     fit_s0 [[j]] <- fit_power_law( sort(-1*IMPACT_s0[[j]],  decreasing=TRUE )*1e-9 )
#     fit_s1a[[j]] <- fit_power_law( sort(-1*IMPACT_s1a[[j]], decreasing=TRUE )*1e-9 )
#     fit_s1b[[j]] <- fit_power_law( sort(-1*IMPACT_s1b[[j]], decreasing=TRUE )*1e-9 )
#     fit_s1c[[j]] <- fit_power_law( sort(-1*IMPACT_s1c[[j]], decreasing=TRUE )*1e-9 )
#       
#     # attach to list
#     CC_list[[j]] <- CC_s1
# 
#   }
# 
#   save( CC_list, IMPACT_s0, IMPACT_s1a, IMPACT_s1b, IMPACT_s1c, fit_s0, fit_s1a, fit_s1b, fit_s1c, file=outfile )
# 
# } else {
# 
#   load(outfile)
# 
# }
# 
# ##------------------------------------------------------------
# ## Get size of events (impact), dataframe 
# ##------------------------------------------------------------
# ## Size difference of the 10 largest extremes
# ## sort by descending s1c impacts
# print("constructing impacts dataframe...")
# list_impacts <- list()
# for (icont in 1:nconts){
#   df_impacts <- tibble( s0 = IMPACT_s0[[icont]], s1a = IMPACT_s1a[[icont]], s1b = IMPACT_s1b[[icont]], s1c = IMPACT_s1c[[icont]] ) %>%
#     filter( s0<0 & s1a<0 & s1b<0 & s1c<0 ) %>%
#     mutate( rank = rank( s1b ),
#             ampl_s1a = s1a / s0, 
#             ampl_s1b = s1b / s0, 
#             ampl_s1c = s1c / s0,
#             icont = icont,
#             event_no = 1:nrow(.)
#           ) %>%
#     mutate( ampl_mean = rowMeans( dplyr::select(., starts_with("ampl")) ) )
#   list_impacts[[icont]] <- df_impacts
# }
# 
# ## combine list into single dataframe (tibble) and get rank by size across all continents (impact in s1b)
# df_impacts <- bind_rows( list_impacts ) %>% 
#               mutate( rank_global = rank( s1b ) ) %>%
#               arrange( s1b )
#               

##------------------------------------------------------------
## Locate events in time and space
##------------------------------------------------------------
## integrate over ten largest events globally
## get dates from NetCDF file
print("plot events...")
nc <- nc_open( "~/data/pmodel_fortran_output/v2/gpp_pmodel_s0_MON_ANOM.nc" )
nc_close(nc)
date <- ymd( "2001-01-01" ) + days( floor(nc$dim$time$vals) )
nevents <- 141
nevents_plot <- 8
cols <- colorRampPalette( brewer.pal(8,"Dark2"))(nevents_plot)

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
    icont <- tmp$icont
    event_no <- tmp$event_no
    rank_global <- tmp$rank_global
    
    ## create array that contains anomaly for voxels in this event and NA otherwise
    ## ... based on s0
    LargestE_s0 <- NA * ANOM_s0.a
    LargestE_s0[CC_list[[icont]]$list[[event_no]]] <- ANOM_s0.a[CC_list[[icont]]$list[[event_no]]]
    
    ## ... based on s1b
    LargestE_s1b <- NA * ANOM_s1b.a
    LargestE_s1b[CC_list[[icont]]$list[[event_no]]] <- ANOM_s1b.a[CC_list[[icont]]$list[[event_no]]]
    
    ## integrate over lon-lat, projecting onto time
    event_time_s0    <- apply( LargestE_s0 , 3, sum, na.rm = TRUE ) * 1e-15
    event_time_s1b   <- apply( LargestE_s1b, 3, sum, na.rm = TRUE ) * 1e-15

    addrows <- tibble( date=date, year=year(date), month=month(date), impact_s0=event_time_s0, impact_s1b=event_time_s1b, icont=icont, event_no=event_no, rank_global=rank_global ) %>%
               filter( abs(impact_s0)>0 | abs(impact_s1b)>0 )
    df_event_time <- bind_rows( df_event_time, addrows )

    ## integrate over time, projecting into lon-lat, based on s1b
    arr_event_lonlat[,,irank] <- apply( LargestE_s1b, 1:2, sum, na.rm = TRUE )
    arr_event_lonlat[ arr_event_lonlat==0 ] <- NA

    ## plot map of impact, integrated over time
    yearchar <- date[which(abs(event_time_s0)>0)] %>% year() %>% table() %>% sort() %>% names()
    if (irank <= nevents_plot) plot_map( -arr_event_lonlat[,,irank], col = cols[irank], add=addmap, file=paste0( "fig/map_event_", sprintf( "%02d", irank ), ".pdf" ), text=yearchar[1] )

    # addmap <- TRUE
  }
  
  save( arr_event_lonlat, df_event_time, file=filn_located )
  
}

## Europe 2003 event
filter(df_event_time, irank==114 ) %>% summarise( impact_s0 = sum(impact_s0), impact_s1b=sum(impact_s1b) ) %>% print()

## Russia 2010 event
filter(df_event_time, irank==9 ) %>% summarise( impact_s0 = sum(impact_s0), impact_s1b=sum(impact_s1b) ) %>% print()

## plot time series of impact, integrated over lon-lat
pdf( "fig/events_overtime.pdf", width = 8, height = 4)
par( mar=c(5,4,1,1), xaxs="i", yaxs="i", las=1)
with( df_event_time, plot(  date, -s1b_e01, type="n", ylim=c(0,max(-s1b_e01, -s1b_e02, -s1b_e03, -s1b_e04, -s1b_e05, -s1b_e06, -s1b_e07, -s1b_e08)), lwd=2, ylab="Impact (PgC / month)" ) )

## s1b
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e01, rep(0,nrow(df_event_time))), col=add_alpha(cols[1], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e02, rep(0,nrow(df_event_time))), col=add_alpha(cols[2], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e03, rep(0,nrow(df_event_time))), col=add_alpha(cols[3], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e04, rep(0,nrow(df_event_time))), col=add_alpha(cols[4], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e05, rep(0,nrow(df_event_time))), col=add_alpha(cols[5], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e06, rep(0,nrow(df_event_time))), col=add_alpha(cols[6], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e07, rep(0,nrow(df_event_time))), col=add_alpha(cols[7], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s1b_e08, rep(0,nrow(df_event_time))), col=add_alpha(cols[8], 0.3), border = NA ) )

## s0
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e01, rep(0,nrow(df_event_time))), col=add_alpha(cols[1], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e02, rep(0,nrow(df_event_time))), col=add_alpha(cols[2], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e03, rep(0,nrow(df_event_time))), col=add_alpha(cols[3], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e04, rep(0,nrow(df_event_time))), col=add_alpha(cols[4], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e05, rep(0,nrow(df_event_time))), col=add_alpha(cols[5], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e06, rep(0,nrow(df_event_time))), col=add_alpha(cols[6], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e07, rep(0,nrow(df_event_time))), col=add_alpha(cols[7], 0.3), border = NA ) )
with( df_event_time, polygon( c(date, rev(date)), c(-s0_e08, rep(0,nrow(df_event_time))), col=add_alpha(cols[8], 0.3), border = NA ) )

legend("topright", as.character(1:9), lty=1, col=cols, bty = "n")
dev.off()

##------------------------------------------------------------
## Plot
##------------------------------------------------------------
magn <- 1
ncols <- 2
nrows <- 1
widths <- rep(6*magn,ncols)
heights <- rep(5*magn,nrows)
order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

pdf( "fig/extremes_global.pdf", width=sum(widths), height = sum(heights) )

  panel <- layout(
          order,
          widths=widths,
          heights=heights,
          TRUE
          )
  # layout.show( panel )

  ## Plot PDF of x>X
  par( xaxs="r", yaxs="r", las=1, mgp=c(3.5,1,0), mar=c(4.5, 4.5, 2, 1))
  n <- nrow(df_impacts)
  plot(   sort( -df_impacts$s1b[1:n], decreasing=TRUE )*1e-9, (1:n)/sum(1:n), log="xy", ylab="p(x)", xlab="Impact (PgC)", pch=16, col = rgb(1,0,0,1), xlim=c(2e4, 6.5e5) )
  polygon( c( sort( -df_impacts$s1a[1:n], decreasing=TRUE )*1e-9, rev(sort( -df_impacts$s1c[1:n], decreasing=TRUE )*1e-9) ), c( (1:n)/sum(1:n), rev((1:n)/sum(1:n)) ), border = NA, col = rgb(1,0,0,0.3) )
  points( sort( -df_impacts$s0[1:n],  decreasing=TRUE )*1e-9, (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
  legend( "bottomleft", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1 )
  mtext( "a", font=2, adj = 0, line = 0.5, cex = 1.2 )
  
  ## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
  par(las=1)
  myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="ratio s1b/s0", col="grey70", names=cont, ylim=c(0.5,2.5) )
  abline( h=1, lty=3 )
  points(1:5, mean_cont$ampl_s1b, pch=16, col="red")
  mtext( "b", font=2, adj = 0, line = 0.5, cex = 1.2 )
  
dev.off()

cont <- c("NA", "SA", "EA", "AF", "AU")
continent <- c("North America", "South America", "Eurasia", "Africa", "Australia")

print("amplification of the largest extremes, globally:")
mean_glob <- df_impacts %>% summarise( ampl_s1a=mean(ampl_s1a), ampl_s1b=mean(ampl_s1b), ampl_s1c=mean(ampl_s1c) )
median_glob <- df_impacts %>% summarise( ampl_s1a=median(ampl_s1a), ampl_s1b=median(ampl_s1b), ampl_s1c=median(ampl_s1c) )
print( mean_glob )
print( median_glob )

mean_cont <- df_impacts %>% group_by(icont) %>% summarise( ampl_s1a=mean(ampl_s1a), ampl_s1b=mean(ampl_s1b), ampl_s1c=mean(ampl_s1c) ) %>% bind_cols( cont=cont )
median_cont <- df_impacts %>% group_by(icont) %>% summarise( ampl_s1a=median(ampl_s1a), ampl_s1b=median(ampl_s1b), ampl_s1c=median(ampl_s1c) ) %>% bind_cols( cont=cont )
print( mean_cont )
print( median_cont )

## Plot difference s1 - s0, absolute
pdf("fig/extremes_s1-s0_abs_glob.pdf",width=10)
  par(las=1)
  with( df_impacts, plot( -1e-15*(s1b-s0), ylab="s1b-s0 (PgC)", xlab="event rank", pch=16, col = rgb(0,0,0,1)) )
  abline( h=0, lty=3 )
dev.off()

## Plot difference s1 - s0, relative
pdf("fig/extremes_s1-s0_rel_glob.pdf",width=10)
  par(las=1)
  with( df_impacts, plot( (s1b/s0), ylab="ratio s1b/s0", xlab="event rank", pch=16, col = rgb(0,0,0,1)) )
  abline( h=1, lty=3 )
dev.off()

## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
pdf( "fig/extremes_s1-s0_rel_glob.pdf", width=7 )
  par(las=1)
  myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="ratio s1b/s0", col="grey70", names=cont, ylim=c(0.5,2.5) )
  abline( h=1, lty=3 )
  points(1:5, mean_cont$ampl_s1b, pch=16, col="red")
dev.off()


# xxxx

# ## Plot PDF of x>X
# cont <- c("NA", "SA", "EU", "AF", "RU", "AU")
# continent <- c("North America", "South America", "Europe", "Africa", "Russia", "Australia")
# pdf( "fig/extremes.pdf", width=10, height = 7 )
# par( mfrow=c(2,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
# for (k in 1:nconts) {
# # par( mfrow=c(1,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
# # for (k in c(1,3,5)) {

#   # df <- tibble( prob=(1:n)/sum(1:n), 
#   #   s1a=sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, 
#   #   s1b=sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, 
#   #   s1c=sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9
#   #   )
#   n0 <- length(IMPACT_s0[[k]])
#   n1 <- length(IMPACT_s1c[[k]])
#   n <- min(n0,n1)
#   plot(   sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy",ylab="p(x)",xlab="PgC", pch=16, col = rgb(1,0,0,1), axes=FALSE )

#   points(   sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )
#   points(   sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )

#   # polygon( c( sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, rev(sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9) ), c( (1:n)/sum(1:n), rev((1:n)/sum(1:n)) ), border = NA, col = rgb(1,0,0,0.3) )
#   points( sort( -1 * IMPACT_s0[[k]][1:n],  decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
#   mtext( bquote( alpha["s1b"] == .(format( fit_s1b[[k]]$alpha, digits = 3) ) ), side=1, line=-4, adj=0.1 )
#   mtext( bquote( alpha["s0"]  == .(format( fit_s0[[k]]$alpha,  digits = 3) ) ), side=1, line=-2, adj=0.1 )
#   mtext( continent[k], line=1, font=2, adj=0 )
#   axis(1, mgp=c(3,1,0) )
#   axis(2, mgp=c(3.5,1,0))
#   box()
# }
# legend( "left", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1.5 )
# dev.off()


# print("amplification of the largest extremes")
# for (k in c(1,3,5,2,4,6)) {
#   ampl <- list_impacts[[k]]$ampl_mean[1:10]
#   print( paste( continent[k], "mean", mean(ampl) ) )
# }
# 
# ## Plot difference s1 - s0, absolute
# pdf("fig/extremes_s1-s0_abs.pdf",width=10)
# par( mfrow=c(2,3), las=1 )
# for (k in c(1,3,5,2,4,6)) {
#   plot( -1e-9*( IMPACT_s1b[[k]] - IMPACT_s0[[k]]), ylab="PgC", pch=16, col = rgb(0,0,0,1) )
#   mtext( continent[k], line=1, font=2, adj=0 )
#   abline( h=0 )
# }
# dev.off()
# 
# ## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
# pdf("fig/extremes_s1-s0_rel.pdf",width=10)
# par( mfrow=c(2,3), las=1 )
# for (k in c(1,3,5,2,4,6)) {
#   myboxplot( list_impacts[[k]]$ampl_mean, ylim=c(0,3), ylab="ratio s1b/s0", col="grey70" )
#   mtext( continent[k], line=1, font=2, adj=0 )
#   abline( h=1, lty=3 )
# }
# dev.off()
# 
# # pdf("fig/extremes_s1-s0_rel_hist.pdf",width=10)
# # par( mfrow=c(2,3), las=1 )
# # for (k in c(1,3,5,2,4,6)) {
# #   hist( ( IMPACT_s1b[[k]] / IMPACT_s0[[k]] ), main="", xlim=c(0,6), breaks=30 )
# #   mtext( continent[k], line=1, font=2, adj=0 )
# # }
# # dev.off()
# 
# 
# # CC <- CC_s1
# # sorted <- sort(unique(as.vector(CC$size)),decreasing=TRUE)
# # ids <- which(CC$size==sorted[1])
# # tmp <- EXT05_s1c*0
# # tmp[ids] <- 1
# # image.plot(apply(tmp,1:2,sum))
# # plot(apply(tmp,3,sum,na.rm=TRUE))
# # 
# # ANOM_s0.a <- aperm(as.array(ANOM_s0),c(2,1,3))[,360:1,]
# # 
# # str(tmp)
# # # has to be summed over the indices...
# # SIZE <- ANOM_s0.a * (CC$size > 0)
# # sorted.real <- sort(unique(as.vector(SIZE)),decreasing=TRUE)
# # 
# # (alpha_s1[k]-1) * fit1$xmin^(alpha_s1[k]-1)
# # x1 <- seq(fit1$xmin, 2, length.out = 100)
# # y1 <- x^(-alpha_s1[k]) * (alpha_s1[k]-1) / sum(x^(-alpha_s1[k]) * (alpha_s1[k]-1))
# # x0 <- seq(fit0$xmin, 1, length.out = 100)
# # y0 <- x^(-alpha_s0[k]) * (alpha_s0[k]-1) / sum(x^(-alpha_s0[k]) * (alpha_s0[k]-1))
# # plot(x1,y1,log="xy")
# # 
# # # playground
# # sorted_s1 <- sort(-1*impact_s1c,decreasing=TRUE, index.return=TRUE)
# # plot(sorted_s1$x*1e-9 + impact_s0[sorted_s1$ix]*1e-9,ylab="PgC/yr",xlab="rank")
# # 
# # plot(-1e-9*(impact_s1c-impact_s0),ylab="PgC/yr")
# # 
# # hist(sorted_s1$x*1e-9 + impact_s0[sorted_s1$ix]*1e-9)
# # abline(h=0)
