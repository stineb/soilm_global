.libPaths( c( .libPaths(), "/home/bstocker/R/x86_64-pc-linux-gnu-library/3.3") )

syshome <- Sys.getenv( "HOME" )
source( paste( syshome, "/.Rprofile", sep="" ) )

library(dplyr)
library(LSD)

siteinfo <- read.csv( paste( myhome, "sofun/input_fluxnet2015_sofun/siteinfo_fluxnet2015_sofun.csv", sep="") )

## Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( paste( "data/nice_all_agg_lue_obs_evi.Rdata", sep="" ) )       # loads 'nice_agg'

##------------------------------------------------
## Bin data w.r.t. alpha
##------------------------------------------------
binwidth <- 0.1
alphabins <- seq( from=0, to=1, by=binwidth )
nice_agg <- nice_agg %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ), insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) 
nice_to_mte_agg <- nice_to_mte_agg %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ), insoilmbin = cut( as.numeric(soilm_mean), breaks = soilmbins ) ) 


par(las=1)
## bias in P-model versus alpha
boxplot( bias_pmodel ~ inalphabin, data=nice_agg, outline=FALSE, col="grey70" )
abline( h=1, lty=2 )

## bias in MTE versus alpha
boxplot( bias_mte ~ inalphabin, data=nice_to_mte_agg, outline=FALSE, col="grey70" )
abline( h=1, lty=2 )

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