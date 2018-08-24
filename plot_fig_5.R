require(fields) ## for image.plot() function
require(dplyr)
require(ncdf4)
require(lubridate)
require(fields)
require(sp)
require(maptools)
library(poweRlaw) # https://cran.r-project.org/web/packages/poweRlaw/vignettes/b_powerlaw_examples.pdf
library(RColorBrewer)
library(igraph) # for fit_power_law()


myboxplot <- function( ... ){
  boxplot( ..., staplewex=0, whisklty=1, outpch=16, outcex=0.5 )
}


load("data/extremes.Rdata")
load("data/impacts_extremes.Rdata")
load("data/extremes_located.Rdata")

##------------------------------------------------------------
## Plot by continent
##------------------------------------------------------------
  cont <- c("NA", "SA", "EA", "AF", "AU")
  continent <- c("North America", "South America", "Eurasia", "Africa", "Australia")
  
  magn <- 0.7
  ncols <- 3
  nrows <- 2
  widths <- rep(6*magn,ncols)
  heights <- rep(5*magn,nrows)
  order <- matrix( seq(nrows*ncols), nrows, ncols, byrow=TRUE)
  
  fct_powerlaw <- function( x, xmin, alpha ){
    xmin * x ^ alpha
  }

  pdf( "fig/fig_5.pdf", width=sum(widths), height = sum(heights) )

    panel <- layout(
            order,
            widths=widths,
            heights=heights,
            TRUE
            )
    # layout.show( panel )

    ## Plot PDF of x>X
    for (icont in 1:5){

      ## load fit data if available, otherwise get new
      filn_powerlawfit <- paste0( "data/d_cpl_cont", icont, ".Rdata" )

      if (file.exists(filn_powerlawfit)){
        
        load( filn_powerlawfit )
      
      } else {

        ## Power Law fitting using "poweRlaw" package
        ## continuous power-law
        d_cpl_s1b = conpl$new(-list_impacts[[icont]]$s1b*1e-15)
        d_cpl_s0  = conpl$new(-list_impacts[[icont]]$s0*1e-15)
        
        ## infer model parameters and update object
        est_s1b = estimate_xmin(d_cpl_s1b)
        if (cont[icont]=="EA"){
          est_s0 <- est_s1b
        } else {
          est_s0  = estimate_xmin(d_cpl_s0 )
        }
        
        d_cpl_s1b$setXmin(est_s1b)
        d_cpl_s0$setXmin(est_s0)

        save( d_cpl_s1b, d_cpl_s0, file = filn_powerlawfit )

      }

      par( xaxs="r", yaxs="r", las=1, mgp=c(3.5,1,0), mar=c(4.5, 4.5, 2, 1))
      
      out <- plot( d_cpl_s1b, pch=16, col = rgb(1,0,0,0.5), xlim=c(range(c(-list_impacts[[icont]]$s1b*1e-15, -list_impacts[[icont]]$s0*1e-15))), xlab="", ylab="" )
      mtext( "Impact (Pg C)", side = 1, line = 3, cex = 0.8 )
      mtext( "p( Impact > x )", side = 2, line = 3.2, las=0, cex = 0.8 )
      out_lines_s1b <- lines(d_cpl_s1b, col = rgb(1,0,0,1), xlab="Impact (PgC)", ylab="p(Impact>x)")
      par(new=TRUE)
      plot( d_cpl_s0, pch=16, col = rgb(0,0,0,0.5), xlim=c(range(c(-list_impacts[[icont]]$s1b*1e-15, -list_impacts[[icont]]$s0*1e-15))), axes=FALSE, xlab="", ylab="" )
      out_lines_s0 <- lines( d_cpl_s0, col = rgb(0,0,0,1))
      
      if (icont==1) legend( "bottomleft", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1.2 )
      mtext( paste0( letters[icont], ") ", continent[icont]), font=2, adj = 0, line = 0.5, cex = 1 )

      ## back-calculate xmin and alpha from lines outout
      alpha_s0 <- (log(out_lines_s0$y[2]) - log(out_lines_s0$y[1])) / (log(out_lines_s0$x[2]) - log(out_lines_s0$x[1]))
      xmin_s0  <- out_lines_s0$y[1] / out_lines_s0$x[1] ^ alpha_s0

      alpha_s1b <- (log(out_lines_s1b$y[2]) - log(out_lines_s1b$y[1])) / (log(out_lines_s1b$x[2]) - log(out_lines_s1b$x[1]))
      xmin_s1b  <- out_lines_s1b$y[1] / out_lines_s1b$x[1] ^ alpha_s1b
      
      pred_s0 <-  tibble( x = seq( min(out_lines_s0$x), max(out_lines_s0$x), length.out = 100 ) ) %>%
                  # mutate( y = fct_powerlaw( x, d_cpl_s0$xmin, -(d_cpl_s0$pars-1) ) )
                  mutate( y = fct_powerlaw( x, xmin_s0, -(d_cpl_s0$pars-1) ) ) %>%  # this works perfectly, but don't know why
                  mutate( y1 = fct_powerlaw( x, xmin_s1b, -(d_cpl_s1b$pars-1) ) )

      pred_s1b <- tibble( x = seq( min(out_lines_s1b$x), max(out_lines_s1b$x), length.out = 100 ) ) %>%
                  # mutate( y = fct_powerlaw( x, d_cpl_s1b$xmin, -(d_cpl_s1b$pars-1) ) )
                  mutate( y = fct_powerlaw( x, xmin_s1b, -(d_cpl_s1b$pars-1) ) ) %>%  # this works perfectly, but don't know why
                  mutate( y0 = fct_powerlaw( x, xmin_s0, -(d_cpl_s0$pars-1) ) )

      lines( out_lines_s0, lwd=2, col="black" )
      lines( out_lines_s1b, lwd=2, col="red" )

      lines( pred_s0, lwd=2, lty=2 )
      lines( pred_s1b, lwd=2, lty=2, col="red" )
      # lines( pred_s1b$x, pred_s1b$y0, lwd=2, lty=2, col="green" )
      
      ## get mean amplification of probability
      # pred_s1b <- pred_s1b %>% mutate( ampl_prob = y/y0 )
      # out <- exp( mean( log(pred_s1b$ampl_prob) ) )
      
      pred_s0 <- pred_s0 %>% mutate( ampl_prob = y1/y )
      ampl <- mean(pred_s0$ampl_prob) # exp( mean( log(pred_s0$ampl_prob) ) )

      print( paste( cont[icont], " Amplification of probability:", ampl ) )

      # mtext( expression( paste( alpha[s0], "= -", d_cpl_s0$pars ) ), line-1, adj = 1 )
      mtext( bquote( italic(N) == .( format( length(list_impacts[[icont]]$s1b), digits = 3 ) ) ), line = -1.4, adj = 0.95, cex=0.9 )
      mtext( bquote( alpha[s0] == .( format( d_cpl_s0$pars, digits = 3 ) ) ), line = -2.7, adj = 0.95, cex=0.9 )
      mtext( bquote( alpha[s1b] == .( format( d_cpl_s1b$pars, digits = 3 ) ) ), line = -3.9, adj = 0.95, cex=0.9, col="red" )
      mtext( bquote( italic(A) == .( format( ampl, digits = 3 ) ) ), line = -5.1, adj = 0.95, cex=0.9, col="red" )
    
    }
    
    ## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
    par(las=1)
      
    myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="s1b/s0", col="grey70", names=cont, ylim=c(0.2,3) )
    abline( h=1, lty=3 )

    mtext( paste0( letters[icont+1], ")"), font=2, adj = 0, line = 0.5, cex = 1 )
    
  dev.off()
