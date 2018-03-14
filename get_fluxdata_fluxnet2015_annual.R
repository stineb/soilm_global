get_fluxdata_fluxnet2015_annual <- function( sitename, add_swcvars=FALSE ){
  ##--------------------------------------------------------------------
  ## Function returns a dataframe containing all the data of flux-derived
  ## GPP for the station implicitly given by path (argument).
  ## Specific for FLUXNET 2015 data
  ## Returns variables in the following units:
  ## temp: deg C
  ## vpd : Pa
  ## prec: mm d-1
  ## nrad: J m-2 d-1
  ## swin: J m-2 d-1
  ## ppfd: mol m-2 d-1 
  ##--------------------------------------------------------------------
  require(dplyr)
  require(readr)
  require(lubridate)

  # ## xxx debug -------------
  # path = "/Users/benjaminstocker/data/FLUXNET-2015_Tier1/20160128/point-scale_none_1d/original/unpacked/FLX_AR-SLu_FLUXNET2015_FULLSET_DD_2009-2011_1-3.csv"
  # add_swcvars = TRUE
  # ## -----------------------

  myhome <- "/alphadata01/bstocker/"

  print( paste( "getting FLUXNET 2015 data for site", sitename ) )
  dirnam_obs <- paste0( myhome, "data/FLUXNET-2015_Tier1/20160128/point-scale_none_1y/original/unpacked/" )
  allfiles <- list.files( dirnam_obs )
  allfiles <- allfiles[ which( grepl( "FULLSET", allfiles ) ) ]
  allfiles <- allfiles[ which( grepl( "3.csv", allfiles ) ) ]
  filnam_obs <- allfiles[ which( grepl( sitename, allfiles ) ) ]
  path <- paste0( dirnam_obs, filnam_obs ) 

  ## from flux to energy conversion, umol/J (Meek et al., 1984), same as used in SPLASH (see Eq.50 in spash_doc.pdf)
  kfFEC <- 2.04

  ## molar mass of C
  c_molmass <- 12.0107

  ## get data
  adf <-  read_csv( path, na="-9999", col_types = cols() ) %>%
          mutate( year = TIMESTAMP )

  ## convert units. given in umolCO2 m-2 s-1. converted to gC m-2 d-1
  adf <- adf %>% mutate(
                        GPP_NT_VUT_REF     = as.numeric(GPP_NT_VUT_REF)    ,
                        GPP_NT_VUT_USTAR50 = as.numeric(GPP_NT_VUT_USTAR50),
                        GPP_DT_VUT_REF     = as.numeric(GPP_DT_VUT_REF)    ,
                        GPP_DT_VUT_USTAR50 = as.numeric(GPP_DT_VUT_USTAR50),
                        LE_F_MDS           = as.numeric(LE_F_MDS),             ## W m-2 -> J m-2 d-1
                        gpp_obs            = ( GPP_NT_VUT_REF + GPP_DT_VUT_REF ) / 2
                        )

  return( adf )
}