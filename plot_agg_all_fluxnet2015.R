.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

library(dplyr)
library(LSD)

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

##------------------------------------------------
## Select only sites that were in NN FLUXNET 2015 analysis
##------------------------------------------------
# Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( paste( "data/nice_all_agg_lue_obs_evi_L2.Rdata", sep="" ) )       # loads 'nice_agg'
load( paste( "data/nice_all_mte_agg_lue_obs_evi_L2.Rdata", sep="" ) )   # loads 'mte_agg'
load( paste( "data/nice_all_modis_agg_lue_obs_evi_L2.Rdata", sep="" ) ) # loads 'modis_agg'

successcodes <- read.csv( paste( myhome, "sofun/utils_sofun/analysis_sofun/fluxnet2015/successcodes.csv", sep="" ), as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename
nice_agg <- nice_agg %>% filter( mysitename %in% do.sites )
mte_agg  <- mte_agg  %>% filter( mysitename %in% do.sites )
modis_agg<- modis_agg%>% filter( mysitename %in% do.sites )


##------------------------------------------------
## Bin data w.r.t. alpha
##------------------------------------------------
binwidth <- 0.2
alphabins <- seq( from=0, to=1, by=binwidth )
soilmbins <- seq( from=0, to=1, by=binwidth )
nice_agg <- nice_agg %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ), insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) 
mte_agg  <- mte_agg  %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ), insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) 
modis_agg<- modis_agg%>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ), insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) 

## get additional variables
cutoff <- 0.5
nice_agg <- nice_agg %>% mutate( dry = ifelse(alpha<cutoff, TRUE, FALSE) )
mte_agg  <- mte_agg  %>% mutate( dry = ifelse(alpha<cutoff, TRUE, FALSE) )
modis_agg<- modis_agg%>% mutate( dry = ifelse(alpha<cutoff, TRUE, FALSE) )

par(las=1)

boxplot( log(bias_modis )  ~ dry, data=modis_agg, outline=FALSE, col="grey70", ylab="log of bias (mod/obs)", xlab=paste("AET/PET <", cutoff), main="MODIS" , xlim=c(0.5,4.5), at=c(1,3))       
boxplot( log(bias_modis * flue_est_nls)  ~ dry, data=modis_agg, outline=FALSE, at=c(2,4), add=TRUE, col=add_alpha("springgreen", 0.5) )
abline( h=0, lty=3 )

boxplot( log(bias_mte )  ~ dry, data=mte_agg, outline=FALSE, col="grey70", ylab="log of bias (mod/obs)", xlab=paste("AET/PET <", cutoff), main="MTE" , xlim=c(0.5,4.5), at=c(1,3))       
boxplot( log(bias_mte * flue_est_nls)  ~ dry, data=mte_agg, outline=FALSE, at=c(2,4), add=TRUE, col=add_alpha("springgreen", 0.5) )
abline( h=0, lty=3 )

boxplot( log(bias_pmodel )  ~ dry, data=nice_agg, outline=FALSE, col="grey70", ylab="log of bias (mod/obs)", xlab=paste("AET/PET <", cutoff), main="P-model" , xlim=c(0.5,4.5), at=c(1,3))       
boxplot( log(bias_pmodel * flue_est_nls)  ~ dry, data=nice_agg, outline=FALSE, at=c(2,4), add=TRUE, col=add_alpha("springgreen", 0.5) )
abline( h=0, lty=3 )


# ## bias in P-model versus alpha
# boxplot( log(bias_pmodel) ~ inalphabin, data=nice_agg, outline=FALSE, col="grey70", main="P-model" )
# abline( h=0, lty=3 )

# ## bias in MTE versus alpha
# boxplot( log(bias_mte) ~ inalphabin, data=mte_agg, outline=FALSE, col="grey70", main="FLUXCOM MTE" )
# abline( h=0, lty=3 )

# ## bias in MODIS versus alpha
# boxplot( log(bias_modis) ~ inalphabin, data=modis_agg, outline=FALSE, col="grey70", main="MODIS" )
# abline( h=0, lty=3 )

# xlim <- c(0,1.2)
# ylim <- c(0,2.5)
# with( 
#       dplyr::filter( nice_agg, ratio_obs_mod<5 & alpha < 0.999 ),  # necessary to get useful bins with heatscatter()
#       heatscatter( 
#                   flue_est, 
#                   ratio_obs_mod, 
#                   xlab="AET/PET",
#                   ylab="GPP observed / GPP modelled",
#                   xlim=xlim,
#                   ylim=ylim,
#                   main=""
#                 )
#     )

# abline( h=1.0, lwd=0.5, lty=2 )
# abline( v=1.0, lwd=0.5, lty=2 )
# lines( c(-99,99), c(-99,99), col='red' )