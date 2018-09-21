source("../utilities/myboxplot.R")

siteinfo <- read.csv( "siteinfo_fluxnet2015_sofun.csv" )

##------------------------------------------------
## Select only sites that were in NN FLUXNET 2015 analysis
##------------------------------------------------
# Load aggregated data from all sites, created by plot_nn_fVAR_fluxnet2015.R: 
load( "data/nice_all_agg_lue_obs_evi.Rdata" )      # loads 'nice_agg', for open-access reproducability, use reduced dataset 'gpp_daily_fluxnet_stocker18natgeo.csv' available from Zenodo XXX instead
load( "data/nice_all_8d_agg_lue_obs_evi.Rdata" )   # loads 'nice_8d', for open-access reproducability, use reduced dataset 'gpp_8daily_fluxnet_stocker18natgeo.csv' available from Zenodo XXX instead

successcodes <- read.csv( "successcodes.csv" ), as.is = TRUE )
do.sites <- dplyr::filter( successcodes, successcode==1 | successcode==2 )$mysitename

## Use only sites where NN method worked (i.e. that had clear and identifiable soil moisture limitation)
nice_agg     <- nice_agg     %>% filter( mysitename %in% do.sites )
nice_8d_agg  <- nice_8d_agg  %>% filter( mysitename %in% do.sites )

##------------------------------------------------
## Bin data w.r.t. alpha
##------------------------------------------------
binwidth <- 0.2
alphabins <- seq( from=0, to=1, by=binwidth )
soilmbins <- seq( from=0, to=1, by=binwidth )
nice_agg     <- nice_agg     %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ) ) 
nice_8d_agg  <- nice_8d_agg  %>% mutate( inalphabin = cut( as.numeric(alpha), breaks = alphabins ) ) 

##------------------------------------------------
## Get bias
##------------------------------------------------
nice_8d_agg <- nice_8d_agg %>%  mutate( bias_pmodel_diff  = gpp_pmodel  - gpp_obs,
                                        bias_modis_diff   = gpp_modis   - gpp_obs,
                                        bias_vpm_diff     = gpp_vpm     - gpp_obs,
                                        bias_bess_v1_diff = gpp_bess_v1 - gpp_obs,
                                        bias_mte_diff     = gpp_mte     - gpp_obs
                                       )
nice_agg <- nice_agg %>%  mutate( bias_pmodel_diff  = gpp_pmodel  - gpp_obs,
                                  bias_bess_v1_diff = gpp_bess_v1 - gpp_obs
                                 )

##------------------------------------------------
## Bin data w.r.t. alpha
##------------------------------------------------

plot_bias_all <- function( nice_agg, nice_8d_agg, filn=NA ){

  xlim <- c(0.5,5.5)
  ylim <- c(-4,4)

  if (!is.na(filn)) pdf(filn, width = 7, height = 6)
    par(xaxs="i", yaxs="i", mgp=c(2.5,1,0), las=1)

    plot( xlim, ylim, type="n", ylim=ylim, xlab = "AET/PET bin", ylab = expression( paste("modelled / observed GPP (ratio)" ) ), xlim=xlim, axes=FALSE )

    rect( 5:1-0.5, rep(ylim[1], 6), 5:1+0.5, rep(ylim[2], 6), border = NA, col=colorRampPalette( c("wheat3", "white") )( 5 ) )

    ## bias in P-model versus alpha
    bp <- myboxplot( log( bias_bess_v1 ) ~ inalphabin, data=nice_agg, add=TRUE, at=5:1+0.3, col="royalblue3", axes=FALSE, boxwex=0.2, outline=FALSE )
    # myboxplot( bias_bess_v1_diff ~ inalphabin, data=nice_agg, add=TRUE, at=5:1+0.3, col="royalblue3", axes=FALSE, boxwex=0.2, outline=FALSE)
    
    text( x = 5:1, y = 3.5, labels = paste( "N =", bp$n ), cex=0.8 )

    ## bias in MODIS versus alpha
    myboxplot( log(bias_vpm) ~ inalphabin, data=nice_8d_agg, add=TRUE, at=5:1+0.1, col="springgreen3", axes=FALSE, boxwex=0.2, outline=FALSE )
    # myboxplot( bias_vpm_diff ~ inalphabin, data=nice_8d_agg, add=TRUE, at=5:1+0.1, col="springgreen3", axes=FALSE, boxwex=0.2, outline=FALSE )
    
    ## bias in VPM versus alpha
    myboxplot( log(bias_modis) ~ inalphabin, data=nice_8d_agg, add=TRUE, at=5:1-0.1, col="orchid", axes=FALSE, boxwex=0.2, outline=FALSE )
    # myboxplot( bias_modis_diff ~ inalphabin, data=nice_8d_agg, add=TRUE, at=5:1-0.1, col="orchid", axes=FALSE, boxwex=0.2, outline=FALSE )
    
    ## bias in BESS v1 versus alpha
    myboxplot( log( bias_pmodel ) ~ inalphabin, data=nice_agg, at=5:1-0.3, add=TRUE, col="tomato", boxwex=0.2, outline=FALSE )
    # myboxplot( bias_pmodel_diff ~ inalphabin, data=nice_agg, at=5:1-0.3, add=TRUE, col="tomato", boxwex=0.2, outline=FALSE )
    
    # ## MTE
    # myboxplot( bias_mte_diff ~ infbin, data = df_dday_8d_agg, add=TRUE, at=5:1-0.4, col="tomato", axes=FALSE, boxwex=0.2 )

    abline( h=0, lty=3 )
    legend("bottomright", c("P-model", "MODIS", "VPM", "BESS"), fill=c("tomato", "orchid", "springgreen3", "royalblue3"), bty="n")

  if (!is.na(filn)) dev.off()
}



