# # make anomalies
# require(raster)
# require(neuroim)
# require(igraph)
# require(fields) ## for image.plot() function
# require(dplyr)
# require(ncdf4)
# require(lubridate)
# 
# myboxplot <- function( ... ){
#   boxplot( ..., staplewex=0, whisklty=1, outline=FALSE )
# }
# 
# ##-----------------------
# overwrite <- FALSE
# ##-----------------------
# 
# outfile <- " "
# 
# if (!file.exists(outfile)||overwrite){
# 
#   ## Define continents from SREX regions
#   # SREX <- raster("/net/exo/landclim/data/dataset/SREX-Region-Masks/20120709/0.5deg_lat-lon_time-invariant/processed/netcdf/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
#   SREX <- raster("~/data/landmasks/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
#   subsets <- list()
#   subsets[[1]] <- 1:6
#   subsets[[2]] <- 7:10
#   subsets[[3]] <- 11:13
#   subsets[[4]] <- 14:17
#   subsets[[5]] <- 18:23
#   subsets[[6]] <- 24:26
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
#   for (j in 1:6) {
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
#     q_all <- quantile( as.vector(ANOM_s1c), c(0.05,0.1,0.9,0.95), na.rm = TRUE )
#     
#     # same quantiles (computed over the combined dataset)
#     EXT05_s1c <- ANOM_s1c < q_all[1]
#     
#     EXT05_s1c <- aperm( as.array(EXT05_s1c), c(2,1,3) )[,latx:1,] == TRUE
#     CC_s1 <- connComp3D( EXT05_s1c )
#     
#     ## keep only events gt 200
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
#     # pdf("fig/extremes.pdf")
#     # plot(   sort( -1*impact_s0,  decreasing=TRUE ) * 1e-9, log="xy", ylab="PgC/yr", xlab="rank" )
#     # points( sort( -1*impact_s1a, decreasing=TRUE ) * 1e-9, col=2 )
#     # points( sort( -1*impact_s1b, decreasing=TRUE ) * 1e-9, col=2 )
#     # points( sort( -1*impact_s1c, decreasing=TRUE ) * 1e-9, col=2 )
#     # legend( "topright", pch=1, col=1:2, legend=c("s0","s1"), bty="n" )
#     # title( paste( "continent", as.character(j) ) )
#     # dev.off()
#     
#     # tmp <- sort( -1*IMPACT_s0[[j]], decreasing=TRUE, index.return = TRUE )
#     # idx_largest[[j]] <- tmp$ix[1]
# 
#     fit_s0 [[j]] <- fit_power_law( tmp$x*1e-9 )
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
# ## Size difference of the 10 largest extremes
# ## sort by descending s1c impacts
# list_impacts <- list()
# for (icont in 1:6){
#   df_impacts <- tibble( s0 = IMPACT_s0[[icont]], s1a = IMPACT_s1a[[icont]], s1b = IMPACT_s1b[[icont]], s1c = IMPACT_s1c[[icont]] ) %>%
#     # arrange( s1c ) %>% 
#     mutate( rank = sort( -1*s1b, decreasing=TRUE, index.return = TRUE )$ix ) %>%
#     mutate( ampl_s1a = s1a / s0, ampl_s1b = s1b / s0, ampl_s1c = s1c / s0 ) %>%
#     mutate( ampl_mean = rowMeans( dplyr::select(., starts_with("ampl")) ) )
#   list_impacts[[icont]] <- df_impacts
# }

## integrate events to get time series
## get dates from NetCDF file
nc <- nc_open( "~/data/pmodel_fortran_output/v2/gpp_pmodel_s0_MON_ANOM.nc" )
nc_close(nc)
date <- ymd( "2001-01-01" ) + days( floor(nc$dim$time$vals) )
cols <- colorRampPalette( brewer.pal(5,"Set1"))(5)

list_event_time   <- list()
list_event_lonlat <- list()

for (icont in 2:2) {

  df_event_time <- tibble( date=date, year=year(date), month=month(date) )
  arr_event_lonlat <- array( NA, dim = c(720,360,5) )
  addmap <- FALSE
  
  for (irank in 1:5){
    
    idx_rank <- which( irank == list_impacts[[icont]]$rank )
    LargestE <- NA * ANOM_s0.a
    LargestE[CC_list[[icont]]$list[[idx_rank]]] <- ANOM_s0.a[CC_list[[icont]]$list[[idx_rank]]]

    ## integrate over lon-lat, projecting onto time
    event_time   <- apply( LargestE, 3, sum, na.rm = TRUE )
    df_tmp <- tibble( x=event_time ) %>% setNames( paste0("e", sprintf( "%02d", idx_rank ) ) )
    df_event_time <- bind_cols( df_event_time, df_tmp )

    ## integrate over time, projecting into lon-lat
    arr_event_lonlat[,,irank] <- apply( LargestE, 1:2, sum, na.rm = TRUE )
    arr_event_lonlat[ arr_event_lonlat==0 ] <- NA

    ## plot map of impact, integrated over time
    ylim <- c(-60,85)
    lat.labels <- seq(-90, 90, 30)
    lat.short  <- seq(-90, 90, 10)
    lon.labels <- seq(-180, 180, 60)
    lon.short  <- seq(-180, 180, 10)
    a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
    b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )
    out.mycolorbar <- mycolorbar( c( "white", cols[irank] ), lev = c( 0, quantile( arr_event_lonlat[,,irank], 0.01, na.rm=TRUE ), 10 ), plot=FALSE ) 
    
    par( mar=c(3,3,3,1),xaxs="i", yaxs="i",las=1)
    image(
      lon, lat, 
      -arr_event_lonlat[,,irank],
      ylim=c(-60,85), 
      yaxt="n", xaxt="n",
      col=colorRampPalette( c( "white", cols[irank] ) )(10),
      xlab="", ylab="",
      add=addmap
    )
    map( add=TRUE, interior=FALSE, resolution=0, lwd=0.5 )
    axis( 2, at=lat.labels, lab=do.call(expression,a), cex.axis=0.7, lwd=1.5 )
    axis( 2, at=lat.short, lab=F, lwd=1, tck=-0.01 )
    axis( 4, at=lat.labels, lab=F, lwd=1.5 )
    axis( 4, at=lat.short, lab=F, lwd=1, tck=-0.01 )
    axis( 1, at=lon.labels, lab=do.call(expression,b), cex.axis=0.7, lwd=1.5 )
    axis( 1, at=lon.short, lab=F, lwd=1, tck=-0.01 )
    axis( 3, at=lon.labels, lab=F, lwd=1.5 )
    axis( 3, at=lon.short, lab=F, lwd=1, tck=-0.01 )
    addmap <- TRUE
  }

  ## plot time series of impact, integrated over lon-lat
  with( df_event_time, plot(  date, e01, col=cols[1], type="l") )
  with( df_event_time, lines( date, e02, col=cols[2] ) )
  with( df_event_time, lines( date, e03, col=cols[3] ) )
  with( df_event_time, lines( date, e04, col=cols[4] ) )
  with( df_event_time, lines( date, e05, col=cols[5] ) )
  
  list_event_time[[icont]] <- df_event_time

}


# 
# ## Plot PDF of x>X
# cont <- c("NA", "SA", "EU", "AF", "RU", "AU")
# continent <- c("North America", "South America", "Europe", "Africa", "Russia", "Australia")
# pdf( "fig/extremes.pdf", width=10, height = 3.5 )
# # par( mfrow=c(2,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
# # for (k in c(1,3,5,2,4,6)) {
# par( mfrow=c(1,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
# for (k in c(1,3,5)) {
# 
#   # df <- tibble( prob=(1:n)/sum(1:n), 
#   #   s1a=sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, 
#   #   s1b=sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, 
#   #   s1c=sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9
#   #   )
# 
#   n0 <- length(IMPACT_s0[[k]])
#   n1 <- length(IMPACT_s1c[[k]])
#   n <- min(n0,n1)
#   plot(   sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy",ylab="p(x)",xlab="PgC", pch=16, col = rgb(1,0,0,1), axes=FALSE )
# 
#   points(   sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )
#   points(   sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )
# 
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
# 
# 
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
