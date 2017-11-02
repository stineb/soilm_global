remove_outliers_fXX <- function( vec, coef=1.5 ) {
  ## use the command boxplot.stats()$out which use the Tukey’s method to identify the outliers ranged above and below the <coef`>*IQR.
  outlier <- boxplot.stats( vec, coef=coef )$out
  outlier <- outlier[ which( outlier>1.0 | outlier<0.0 ) ]
  vec[ which( is.element( vec, outlier ) ) ] <- NA
  return( vec )
}

remove_outliers <- function( vec, coef=1.5 ) {
  ## use the command boxplot.stats()$out which use the Tukey’s method to identify the outliers ranged above and below the <coef`>*IQR.
  outlier <- boxplot.stats( vec, coef=coef )$out
  vec[ which( is.element( vec, outlier ) ) ] <- NA
  return( vec )
}