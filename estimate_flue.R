library(dplyr)
library(caret)
library(nnet)
library(randomForest)
library(tibble)
library(readr)

source("analyse_modobs.R")

## load dataframe with data aligned by drought event
load("data/data_aligned_agg.Rdata") ## loads 'df_dday_agg'
df <- as_tibble(df_dday_agg)

# ## load dataframe with all data in there
# load("data/nice_nn_agg_lue_obs_evi.Rdata")
# df <- as_tibble(nice_agg)

## add mean alpha (AET/PET) to dataframe
load( "data/alpha_fluxnet2015.Rdata" )  # loads 'df_alpha'
df <- df %>% left_join( rename( df_alpha, meanalpha=alpha ), by="mysitename" )

## add mean aridity index (P/PET) to dataframe
load( "data/ai_fluxnet2015.Rdata" )  # loads 'df_ai'
df <- df %>% left_join( select( df_ai, mysitename, ai ), by="mysitename" )

## add water table depth to dataframe
load( "data/wtd_fluxnet2015.Rdata" )  # loads 'df_wtd'
df <- df %>% left_join( df_wtd, by="mysitename" )

## add vegetation type ('classid') to dataframe
siteinfo <- read_csv( "data/siteinfo_fluxnet2015_sofun.csv" )
df <- df %>% left_join( select( siteinfo, mysitename, classid ), by="mysitename" ) %>%
             mutate( GRA = ifelse( classid=="GRA", TRUE, FALSE ),
                     SAV = ifelse( classid=="SAV", TRUE, FALSE ),
                     ENF = ifelse( classid=="ENF", TRUE, FALSE ),
                     WET = ifelse( classid=="WET", TRUE, FALSE ),
                     WSA = ifelse( classid=="WSA", TRUE, FALSE ),
                     EBF = ifelse( classid=="EBF", TRUE, FALSE ),
                     DBF = ifelse( classid=="DBF", TRUE, FALSE ),
                     CSH = ifelse( classid=="CSH", TRUE, FALSE )
              ) %>%
             mutate(  classid = as.factor(classid), 
                      GRA = as.factor(GRA),
                      SAV = as.factor(SAV),
                      ENF = as.factor(ENF),
                      WET = as.factor(WET),
                      WSA = as.factor(WSA),
                      EBF = as.factor(EBF),
                      DBF = as.factor(DBF),
                      CSH = as.factor(CSH)
                      )

## clean dataframe (drop rows where at least 1 of columns is NA)
df <- df %>% select( mysitename, date, fvar, soilm_splash, classid, wtd, ai, fpar, GRA, SAV, ENF, WET, WSA, EBF, DBF, CSH )

na_idxs <- apply( df, 2, FUN=function (x) which(is.na(x)) ) %>% unlist() %>% unique()
df <- df[-na_idxs,]
print(  apply( df, 2, FUN=function (x) sum(is.na(x)) ) )

set.seed(1)

## training Sample with 80% of the data
train = sample( 1:nrow(df), round(0.8*nrow(df)) )

##------------------------------------------------------------------
## Random Forest
##------------------------------------------------------------------
print("random forest...")
rf = randomForest( 
  # fvar ~ soilm_splash + ai + fpar , # + wtd + GRA
  # fvar ~ soilm_splash + GRA , # + wtd 
  fvar ~ soilm_splash + ai , # + wtd + GRA
  data = df,
  importance = TRUE,
  subset = train 
  )

analyse_modobs( df$fvar, predict( rf, df ) )
varImpPlot(rf)
save( rf, file = "data/rf.Rdata")


##------------------------------------------------------------------
## Neural Network
##------------------------------------------------------------------
print("neural network...")
nn_nnet <- nnet( 
  fvar ~ soilm_splash + classid + ai + wtd + fpar, 
  data = df, 
  linout = TRUE, 
  subset = train,
  size = 20 
  )

par(mfrow=c(1,1))
analyse_modobs( df[train,"fvar"]$fvar, nn_nnet$fitted.values[,1] )
varImp( nn_nnet )
save( nn_nnet, file = "data/nn_nnet.Rdata")

##------------------------------------------------------------------
## Caret
##------------------------------------------------------------------
## scale data to within 0 and 1
preprocessParams <- preProcess( df, method=c("center", "scale") )

## define parameters for training
traincotrlParams <- trainControl( method="repeatedcv", 
                                  number=3, repeats=3, verboseIter=FALSE, p=0.75 
                                )

## sample number of nodes for best performance
tune_grid <- expand.grid( .decay = c(0.1), .size = 16 )  ## found best size to be 16

## train
print("caret nnet...")
nn_caret_classid <- train(
  fvar ~ soilm_splash + classid + ai + wtd + fpar,
  data      = df,
  method    = "nnet",
  tuneGrid  = tune_grid,
  trControl = traincotrlParams,
  trace     = FALSE
)

analyse_modobs( df[,"fvar"]$fvar, as.vector(predict( nn_caret_classid, df )) )
varImp( nn_caret_classid )
save( nn_caret_classid, file = "data/nn_caret_classid.Rdata" )

## one sees here that within the 'classid' predictor, it's the highest importance is whether or not it is a GRAland ("GRA")
## therefore, train again, now with 'GRA' as a predictor, instead of 'classid'
## train
print("caret nnet (GRA)...")
nn_caret_GRA <- train(
  fvar ~ soilm_splash + GRA + ai + fpar,
  data      = df,
  method    = "nnet",
  tuneGrid  = tune_grid,
  trControl = traincotrlParams,
  trace     = FALSE
)

analyse_modobs( df[,"fvar"]$fvar, as.vector(predict( nn_caret_GRA, df )) )
varImp( nn_caret_GRA )
save( nn_caret_GRA, file = "data/nn_caret_GRA.Rdata" )

# ## train
# print("caret rf....")
# rf_caret <- train(
#   fvar ~ soilm_splash + GRA + ai + fpar,
#   data      = df,
#   method    = "rf",
#   trControl = traincotrlParams
# )

# analyse_modobs( df[,"fvar"]$fvar, as.vector(predict( rf_caret, df )) )
# varImp( rf_caret )
# save( rf_caret, file = "data/rf_caret.Rdata" )

##------------------------------------------------------------------
## Greedy search of priority list of predictors
##------------------------------------------------------------------
greedy_pred <- function( df, target, preds ){

  ##///////////////////////////////////////////////////////////////////
  ## This function consecutively adds predictors to NN model finding 
  ## the order of best performance.
  ##-------------------------------------------------------------------
  set.seed(1)

  ## define preprocessing
  preprocessParams <- preProcess( df, method=c("center", "scale") )

  ## define parameters for training
  traincotrlParams <- trainControl( method="repeatedcv", 
                                    number=3, repeats=3, verboseIter=FALSE, p=0.75 
                                  )

  ## sample number of nodes for best performance (don't actually expand anything here)
  tune_grid <- expand.grid( .decay = c(0.1), .size = 16 )  ## found best size to be 16

  use_preds <- c()
  prio <- data.frame()

  for (npreds in 1:length(preds)){

    print("----------------------------------------")
    print( paste( "## Number of predictors:", npreds ) )

    prio_npred <- data.frame()

    if (npreds==1){
      eval_preds <- preds
    } else {
      eval_preds <- preds[ !(preds %in% use_preds) ]
    }
    
    for (ipred in eval_preds){

      use_preds_tmp <- c( use_preds, ipred )
      forml <- paste( target, "~", paste( use_preds_tmp, collapse=" + " ) )
      print( forml )

      ## train
      nn <- train(
        as.formula(forml),
        data      = df,
        method    = "nnet",
        tuneGrid  = tune_grid,
        trControl = traincotrlParams,
        trace     = FALSE
      )
      
      ## get performance statistics of model: modelled (NN) vs. observed
      out <- analyse_modobs( df[[ target ]], as.vector( predict( nn, df ) ), do.plot=FALSE )
      
      ## construct "priority array" = array holding all the performance statistics
      tmp <- data.frame( nse=out$nse, rsq=out$rsq, rmse=out$rmse )
      row.names(tmp) <- ipred
      prio_npred <- rbind( prio_npred, tmp )

    }

    ## determine which configuration was best
    ## Sort prio predictor perfomance data frame by R-squared
    prio_npred <- prio_npred[ order(-prio_npred$rsq), , drop = FALSE ]
    print( prio_npred )

    ## retain only best-performing predictor
    use_preds <- c( use_preds, rownames(prio_npred[1,]) )
    prio <- rbind( prio, prio_npred[1,] )

  }

  out <- list( prio=prio, nn=nn )

  print(out)
  return( out )

}

out <- greedy_pred( df=df, target = "fvar", preds =  c( "soilm_splash", "ai", "wtd", "fpar", "GRA", "SAV" ) ) # , "ENF", "WET", "WSA", "EBF", "DBF", "CSH"
save( out, file="data/prio_preds.Rdata" )

#                      nse       rsq       rmse
# soilm_splash -0.08188449 0.4851471 0.15730475
# classid       0.58598194 0.7139707 0.11728472
# ai            0.71382006 0.7843892 0.10187333
# fpar          0.74363376 0.8047961 0.09699041
# wtd           0.74852281 0.8091091 0.09595733


