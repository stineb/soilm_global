binned_bias_sse <- function( df, par ){

  require(dplyr)
  require(tidyr)

  source("stress_quad_1sided.R")

  ## calculate stress as a function of soil moisture, and meanalpha given apar (par[1]) and bpar (par[2])
  df <- df %>% mutate( flue_est = stress_quad_1sided_alpha( soilm_mean, meanalpha, 0.9, par[1], par[2] ) ) %>%

               ## correct modelled, i.e. modify bias accordingly
               mutate( bias_pmodel_corr = bias_pmodel * flue_est ) 

  ## bin data vs. fLUE (fvar)
  nbins <- 10
  bins <- seq( 0.0, 1.0, 1.0/nbins )
  xvals <- bins[1:(length(bins)-1)] + (bins[2]-bins[1])/2
  df <- tibble( fvar=df$fvar, bias_pmodel_corr=df$bias_pmodel_corr )
  df <- df %>% mutate( inbenibin = cut( fvar , breaks = bins ) ) %>% group_by( inbenibin )
  tmp <- df %>% summarise( median=median( bias_pmodel_corr, na.rm=TRUE ) ) %>% complete( inbenibin, fill = list( median = NA ) ) %>% dplyr::select( median )
  bias_pmodel_corr <- unlist(tmp)[1:nbins]

  ## calculate sum of square errors of median within bins
  sse <- sum( log(bias_pmodel_corr) )

  ## return sse because that's what needs to be minimised
  return(sse)

}

out <- optim( par=c(0,1), binned_bias_sse, df=filter(ddf, bias_pmodel < 10)  )

## plot
## calculate stress as a function of soil moisture, and meanalpha given apar (par[1]) and bpar (par[2])
ddf <- ddf %>% mutate( flue_est = stress_quad_1sided_alpha( soilm_mean, meanalpha, 0.9, out$par[1], out$par[2] ) ) %>%
  
  ## correct modelled, i.e. modify bias accordingly
  mutate( bias_pmodel_corr = bias_pmodel * flue_est ) 

## bin data vs. fLUE (fvar)
nbins <- 5
bins <- seq( 0.0, 1.0, 1.0/nbins )
xvals <- bins[1:(length(bins)-1)] + (bins[2]-bins[1])/2
tmp <- tibble( fvar=ddf$fvar, bias_pmodel_corr=ddf$bias_pmodel_corr )
tmp <- tmp %>% mutate( inbenibin = cut( fvar , breaks = bins ) ) %>% group_by( inbenibin )
boxplot( log(bias_pmodel_corr) ~ inbenibin, data=tmp, outline=FALSE )
abline( h=0, lty=3 )
