# make anomalies
library(raster)
library(neuroim)
library(igraph)

# rasterOptions( tmpdir="/net/exo/landclim/zjakob/tmp", tmptime = 240 )
# fdir <- "/net/firebolt/data/zjakob/Stocker/"
# setwd(fdir)

# gpp_s0 <- brick("s0_gpp_all_detrended.m.nc")
# MSC_s0 <- calc(gpp_s0, function(x) {if (is.na(x[1])) return(rep(NA,12)) else return(rowMeans(matrix(x,12,30)))})

# gpp_s1 <- brick("s1_gpp_all_detrended.m.nc")
# MSC_s1 <- calc(gpp_s1, function(x) {if (is.na(x[1])) return(rep(NA,12)) else return(rowMeans(matrix(x,12,30)))})

# gpp_s1.tmp <- brick("s1_gpp_all.m.nc")

# # compute anomalies
# ANOM_s0 <- gpp_s0 - MSC_s0
# ANOM_s1 <- gpp_s1 - MSC_s1

# # save anomalies
# # writeRaster(ANOM_s0, "/net/firebolt/data/zjakob/Stocker/s0_gpp_anomalies.m.nc", format="CDF", varname="GPP anomalies", overwrite=TRUE)
# # writeRaster(ANOM_s1, "/net/firebolt/data/zjakob/Stocker/s1_gpp_anomalies.m.nc", format="CDF", varname="GPP anomalies", overwrite=TRUE)

##########################################################################################################################
# start from here
# SREX <- raster("/net/exo/landclim/data/dataset/SREX-Region-Masks/20120709/0.5deg_lat-lon_time-invariant/processed/netcdf/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
SREX <- raster("~/data/landmasks/srex-region-masks_20120709.srex_mask_SREX_masks_all.05deg.time-invariant.nc")
subsets <- list()
subsets[[1]] <- 1:6
subsets[[2]] <- 7:10
subsets[[3]] <- 11:13
subsets[[4]] <- 14:17
subsets[[5]] <- 18:23
subsets[[6]] <- 24:26

alpha_s0 <- c()
alpha_s1 <- c()
IMPACT_s0 <- list()
IMPACT_s1 <- list()

ANOM_s0_all <- brick("~/data/pmodel_fortran_output/gpp_pmodel_s0_DAY_ANOM.nc")
ANOM_s1_all <- brick("~/data/pmodel_fortran_output/gpp_pmodel_s1b_DAY_ANOM.nc")

# ANOM_s0.a <- aperm(as.array(ANOM_s0_all * area(ANOM_s0_all)),c(2,1,3))[,360:1,]
# ANOM_s1.a <- aperm(as.array(ANOM_s1_all * area(ANOM_s1_all)),c(2,1,3))[,360:1,]

latx <- 360
for (j in 1:6) {
  
  print(j)
  
  # only 1 continent
  Na <- SREX %in% subsets[[j]]

  # ANOM_s0 <- mask(ANOM_s0_all, Na, maskvalue=0)
  ANOM_s1 <- mask(ANOM_s1_all, Na, maskvalue=0)
  
  # q_all <- quantile(c(as.vector(ANOM_s1), as.vector(ANOM_s0)), c(0.05,0.1,0.9,0.95), na.rm = TRUE)
  q_all <- quantile( as.vector(ANOM_s1), c(0.05,0.1,0.9,0.95), na.rm = TRUE )
  
  # same quantiles (computed over the combined dataset)
  # only for s1
  EXT05_s1 <- ANOM_s1 < q_all[1]
  # EXT05_s0 <- ANOM_s1 < q_all[1]
  
  # EXT05_s0.a <- aperm(as.array(EXT05_s0),c(2,1,3))[,latx:1,] == TRUE
  # CC_s0 <- connComp3D(EXT05_s0.a)
  EXT05_s1.a <- aperm( as.array(EXT05_s1), c(2,1,3) )[,latx:1,] == TRUE
  CC_s1 <- connComp3D( EXT05_s1.a )
  
  # CC_s0$index[CC_s0$size<200] <- 0
  CC_s1$index[CC_s1$size < 200] <- 0
  
  # make list with indices pointing to equal index
  # CC_s0$list <- sapply(1:max(CC_s0$index), function(x) which(CC_s0$index == x))
  CC_s1$list <- sapply( 1:max(CC_s1$index), function(x) which(CC_s1$index == x) )
  
  # find largest extremes defined by impact
  # impact_s0 <- unlist(lapply(CC_s0$list, function(x) sum(ANOM_s0.a[x])))
  impact_s1 <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s1.a[x])) )
  impact_s0 <- unlist( lapply( CC_s1$list, function(x) sum(ANOM_s0.a[x])) )
  
  
  # pdf("/net/firebolt/data/zjakob/Stocker/extremes.pdf")
  # plot(sort(-1*impact_s0,decreasing=TRUE)*1e-9, log="xy",ylab="PgC/yr",xlab="rank")
  # points(sort(-1*impact_s1,decreasing=TRUE)*1e-9,col=2)
  # legend("topright",pch=1,col=1:2,legend=c("s0","s1"),bty="n")
  # dev.off()
  
  IMPACT_s1[[j]] <- impact_s1
  IMPACT_s0[[j]] <- impact_s0
  
  fit1 <- fit_power_law(sort(-1*impact_s1,decreasing=TRUE)*1e-9)
  fit0 <- fit_power_law(sort(-1*impact_s0,decreasing=TRUE)*1e-9)
  
  alpha_s1[j] <- fit1$alpha
  alpha_s0[j] <- fit0$alpha
  
}
str(impact_s0)


fit0 <- fit_power_law(sort(-1*IMPACT_s0[[k]][1:n],decreasing=TRUE)*1e-9)
fit1 <- fit_power_law(sort(-1*IMPACT_s1[[k]][1:n],decreasing=TRUE)*1e-9)


cont <- c("NA", "SA", "EU", "AF", "RU", "AU")
pdf("/net/firebolt/data/zjakob/Stocker/extremes.pdf",width=10)
pdf("/net/firebolt/data/zjakob/Stocker/extremes_s0_based_on_s1.pdf",width=10)
par(mfrow=c(2,3))
for (k in c(1,3,5,2,4,6)) {
  n0 <- length(IMPACT_s0[[k]])
  n1 <- length(IMPACT_s1[[k]])
  n <- min(n0,n1)
  plot(sort(-1*IMPACT_s1[[k]][1:n],decreasing=TRUE)*1e-9, (1:n)/sum(1:n), log="xy",ylab="p(x)",xlab="PgC/yr", main=cont[k])
  points(sort(-1*IMPACT_s0[[k]][1:n],decreasing=TRUE)*1e-9, (1:n)/sum(1:n), col=2)
  mtext(paste("alpha_s1 =", round(alpha_s1[[k]],2)), side=3, line=-2, adj=0.9)
  mtext(paste("alpha_s0 =", round(alpha_s0[[k]],2)), side=3, line=-4, adj=0.9)
  
}
legend("bottomleft",pch=1,col=1:2,legend=c("s1","s0"),bty="n")
dev.off()

pdf("/net/firebolt/data/zjakob/Stocker/extremes_s1-s0.pdf",width=10)
par(mfrow=c(2,3))
for (k in c(1,3,5,2,4,6)) {
  plot(-1e-9*(IMPACT_s1[[k]]-IMPACT_s0[[k]]),ylab="PgC/yr")
  abline(h=0)
}
dev.off()


abline(h=0)


CC <- CC_s1
sorted <- sort(unique(as.vector(CC$size)),decreasing=TRUE)
ids <- which(CC$size==sorted[1])
tmp <- EXT05_s0.a*0
tmp[ids] <- 1
image.plot(apply(tmp,1:2,sum))
plot(apply(tmp,3,sum,na.rm=TRUE))

ANOM_s0.a <- aperm(as.array(ANOM_s0),c(2,1,3))[,360:1,]

str(tmp)
# has to be summed over the indices...
SIZE <- ANOM_s0.a * (CC$size > 0)
sorted.real <- sort(unique(as.vector(SIZE)),decreasing=TRUE)

(alpha_s1[k]-1) * fit1$xmin^(alpha_s1[k]-1)
x1 <- seq(fit1$xmin, 2, length.out = 100)
y1 <- x^(-alpha_s1[k]) * (alpha_s1[k]-1) / sum(x^(-alpha_s1[k]) * (alpha_s1[k]-1))
x0 <- seq(fit0$xmin, 1, length.out = 100)
y0 <- x^(-alpha_s0[k]) * (alpha_s0[k]-1) / sum(x^(-alpha_s0[k]) * (alpha_s0[k]-1))
plot(x1,y1,log="xy")

# playground
sorted_s1 <- sort(-1*impact_s1,decreasing=TRUE, index.return=TRUE)
plot(sorted_s1$x*1e-9 + impact_s0[sorted_s1$ix]*1e-9,ylab="PgC/yr",xlab="rank")

plot(-1e-9*(impact_s1-impact_s0),ylab="PgC/yr")

hist(sorted_s1$x*1e-9 + impact_s0[sorted_s1$ix]*1e-9)
abline(h=0)
