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

plot_map <- function( out, col, add, file=NA, text="" ){
  ylim <- c(-60,85)
  lon <- seq(-179.75, 179.75, 0.5)
  lat  <- seq(-89.75, 89.75, 0.5)
  lat.labels <- seq(-90, 90, 30)
  lat.short  <- seq(-90, 90, 10)
  lon.labels <- seq(-180, 180, 60)
  lon.short  <- seq(-180, 180, 10)
  a <- sapply( lat.labels, function(x) bquote(.(x)*degree ~ N) )
  b <- sapply( lon.labels, function(x) bquote(.(x)*degree ~ E) )
  
  if (!is.na(file)) pdf( file, width=8, height=6 )
    par( mar=c(3,3,3,1),xaxs="i", yaxs="i",las=1)
    image(
      lon, lat, 
      out,
      ylim=c(-60,85), 
      yaxt="n", xaxt="n",
      col=colorRampPalette( c( add_alpha(col, 0.3), col ) )(10),
      xlab="", ylab="",
      add=add
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
    text( -160, -50, text, col=col )
  if (!is.na(file)) dev.off()  

}


load("data/extremes.Rdata")
load("data/impacts_extremes.Rdata")
load("data/extremes_located.Rdata")


print("amplification of the largest extremes, globally:")
mean_glob <- df_impacts %>% summarise( ampl_s1a=mean(ampl_s1a), ampl_s1b=mean(ampl_s1b), ampl_s1c=mean(ampl_s1c), s0=mean(s0), s1a=mean(s1a), s1b=mean(s1b), s1c=mean(s1c) )
median_glob <- df_impacts %>% summarise( ampl_s1a=median(ampl_s1a), ampl_s1b=median(ampl_s1b), ampl_s1c=median(ampl_s1c), s0=median(s0), s1a=median(s1a), s1b=median(s1b), s1c=median(s1c) )
print( mean_glob )
print( median_glob )

cont <- c("NA", "SA", "EA", "AF", "AU")
mean_cont <- df_impacts %>% group_by(icont) %>% summarise( ampl_s1a=mean(ampl_s1a), ampl_s1b=mean(ampl_s1b), ampl_s1c=mean(ampl_s1c), s0=mean(s0), s1a=mean(s1a), s1b=mean(s1b), s1c=mean(s1c) ) %>% bind_cols( cont=cont )
median_cont <- df_impacts %>% group_by(icont) %>% summarise( ampl_s1a=median(ampl_s1a), ampl_s1b=median(ampl_s1b), ampl_s1c=median(ampl_s1c), s0=median(s0), s1a=median(s1a), s1b=median(s1b), s1c=median(s1c) ) %>% bind_cols( cont=cont )
print( mean_cont )
print( median_cont )

## Europe 2003 event
filter(df_event_time, rank_global==15 ) %>% summarise( impact_s0 = sum(impact_s0), impact_s1b=sum(impact_s1b) ) %>% print()

## Russia 2010 event
filter(df_event_time, rank_global==11 ) %>% summarise( impact_s0 = sum(impact_s0), impact_s1b=sum(impact_s1b) ) %>% print()

## plot time series of impact, integrated over lon-lat
nevents_plot <- 8
cols <- colorRampPalette( brewer.pal(8,"Dark2"))(nevents_plot)
pdf( "fig/events_overtime.pdf", width = 8, height = 4)
par( mar=c(5,4,1,1), xaxs="i", yaxs="i", las=1)
with( df_event_time, plot(  date, -impact_s1b, type="n", ylim=c(0,max(-impact_s1b)), lwd=2, ylab="Impact (PgC/month)", xlim=c( ymd("1982-01-21"), ymd("2019-12-12")) ) )
for (irank in 1:nevents_plot){
  with( filter( df_event_time, rank_global==irank ), polygon( c(date, rev(date)), c(-impact_s1b, rep(0,nrow(filter( df_event_time, rank_global==irank )))), col=add_alpha(cols[irank], 0.6), border = NA ) )
}
legend("topright", as.character(1:nevents_plot), fill=add_alpha(cols, 0.6), bty = "n", border = NA)
dev.off()

##------------------------------------------------------------
## Plot Global
##------------------------------------------------------------
  magn <- 1
  ncols <- 2
  nrows <- 1
  widths <- rep(6*magn,ncols)
  heights <- rep(5*magn,nrows)
  order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

  pdf( "fig/extremes_global.pdf", width=sum(widths), height = sum(heights) )

    panel <- layout(
            order,
            widths=widths,
            heights=heights,
            TRUE
            )
    # layout.show( panel )

    ## Plot PDF of x>X
    par( xaxs="r", yaxs="r", las=1, mgp=c(3.5,1,0), mar=c(4.5, 4.5, 2, 1))
    n <- nrow(df_impacts)
    plot(   sort( -df_impacts$s1b[1:n]*1e-15, decreasing=TRUE ), (1:n)/sum(1:n), log="xy", ylab="p(x)", xlab="Impact (PgC)", pch=16, col = rgb(1,0,0,1) )
    polygon( c( sort( -df_impacts$s1a[1:n]*1e-15, decreasing=TRUE ), rev(sort( -df_impacts$s1c[1:n]*1e-15, decreasing=TRUE )) ), c( (1:n)/sum(1:n), rev((1:n)/sum(1:n)) ), border = NA, col = rgb(1,0,0,0.3) )
    points( sort( -df_impacts$s0[1:n]*1e-15,  decreasing=TRUE ), (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
    legend( "bottomleft", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1 )
    mtext( "a", font=2, adj = 0, line = 0.5, cex = 1.2 )
    
    ## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
    par(las=1)
    myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="ratio s1b/s0", col="grey70", names=cont, ylim=c(0.5,2.5) )
    abline( h=1, lty=3 )
    points(1:5, mean_cont$ampl_s1b, pch=16, col="red")
    mtext( "b", font=2, adj = 0, line = 0.5, cex = 1.2 )
    
  dev.off()


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

  pdf( "fig/extremes_bycont.pdf", width=sum(widths), height = sum(heights) )

    panel <- layout(
            order,
            widths=widths,
            heights=heights,
            TRUE
            )
    # layout.show( panel )

    ## Plot PDF of x>X
    for (icont in 1:5){
      
      par( xaxs="r", yaxs="r", las=1, mgp=c(3.5,1,0), mar=c(4.5, 4.5, 2, 1))

      # # n <- nrow(list_impacts[[icont]])
      # # plot(       sort( -list_impacts[[icont]]$s1b[1:n]*1e-15, decreasing=TRUE ), (1:n)/sum(1:n), log="xy", ylab="p(x)", xlab="Impact (PgC)", pch=16, col = rgb(1,0,0,1), xlim=c(range(c(-list_impacts[[icont]]$s1b[1:n]*1e-15, -list_impacts[[icont]]$s0[1:n]*1e-15))) )
      # # # polygon( c( sort( -list_impacts[[icont]]$s1a[1:n]*1e-15, decreasing=TRUE ), rev(sort( -list_impacts[[icont]]$s1c[1:n]*1e-15, decreasing=TRUE )) ), c( (1:n)/sum(1:n), rev((1:n)/sum(1:n)) ), border = NA, col = rgb(1,0,0,0.3) )
      # # points(     sort( -list_impacts[[icont]]$s0[1:n]*1e-15,  decreasing=TRUE ), (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
      # # legend( "bottomleft", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1 )
      # # mtext( paste0( letters[icont], ") ", continent[icont]), font=2, adj = 0, line = 0.5, cex = 1 )
      # 
      # ## no cutoff
      # n <- length(list_impacts[[icont]]$s1b)
      # plot(       sort( -list_impacts[[icont]]$s1b*1e-15, decreasing=TRUE ), (1:n)/sum(1:n), log="xy", ylab="p(x)", xlab="Impact (PgC)", pch=16, col = rgb(1,0,0,1), xlim=c(range(c(-list_impacts[[icont]]$s1b*1e-15, -list_impacts[[icont]]$s0*1e-15))) )
      # points(     sort( -list_impacts[[icont]]$s0*1e-15,  decreasing=TRUE ), (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
      # legend( "bottomleft", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1 )
      # mtext( paste0( letters[icont], ") ", continent[icont]), font=2, adj = 0, line = 0.5, cex = 1 )
      # 
      # ## Power Law fitting using "igraph" package
      # fit0 <- fit_power_law(-list_impacts[[icont]]$s0*1e-15)
      # fit1 <- fit_power_law(-list_impacts[[icont]]$s1b*1e-15)
      # 
      # pred_s0 <-  tibble( x = seq( fit0$xmin, max(-list_impacts[[icont]]$s0*1e-15), length.out = 100 ) ) %>%
      #             mutate( y = fct_powerlaw( x, fit0$xmin, -fit0$alpha ))
      # pred_s1b <- tibble( x = seq( fit1$xmin, max(-list_impacts[[icont]]$s1b*1e-15), length.out = 100 ) ) %>%
      #             mutate( y = fct_powerlaw( x, fit1$xmin, -fit1$alpha ) )
      # 
      # lines( pred_s0, lwd=2 )
      # lines( pred_s1b, lwd=2, col="red" )
      # 
      # # alpha_s1 <- fit1$alpha
      # # alpha_s0 <- fit0$alpha
      # 
      # # mtext(paste("alpha_s1 =", round(alpha_s1[[k]],2)), side=3, line=-2, adj=0.9)
      # # mtext(paste("alpha_s0 =", round(alpha_s0[[k]],2)), side=3, line=-4, adj=0.9)
      # 
      # (alpha_s1-1) * fit1$xmin^(alpha_s1-1)
      # x1 <- seq(fit1$xmin, 2, length.out = 100)
      # y1 <- x1^(-alpha_s1) * (alpha_s1-1) / sum(x1^(-alpha_s1) * (alpha_s1-1))
      # x0 <- seq(fit0$xmin, 1, length.out = 100)
      # y0 <- x0^(-alpha_s0) * (alpha_s0-1) / sum(x0^(-alpha_s0) * (alpha_s0-1))
      # lines(x1,y1,log="xy")

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
      
      out <- plot( d_cpl_s1b, pch=16, col = rgb(1,0,0,0.5), xlim=c(range(c(-list_impacts[[icont]]$s1b*1e-15, -list_impacts[[icont]]$s0*1e-15))), xlab="", ylab="" )
      mtext( "Impact (PgC)", side = 1, line = 3, cex = 0.8 )
      mtext( "p(Impact>x)", side = 2, line = 3.2, las=0, cex = 0.8 )
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
      
      # (alpha_s1-1) * fit1$xmin^(alpha_s1-1)
      # x1 <- seq(fit1$xmin, 2, length.out = 100)
      # y1 <- x1^(-alpha_s1) * (alpha_s1-1) / sum(x1^(-alpha_s1) * (alpha_s1-1))
      # x0 <- seq(fit0$xmin, 1, length.out = 100)
      # y0 <- x0^(-alpha_s0) * (alpha_s0-1) / sum(x0^(-alpha_s0) * (alpha_s0-1))
      # lines(x1,y1, col="blue")
      # 
      # ## test if data is drawn from a power law distribution
      # bs_cpl = bootstrap_p( d_cpl, no_of_sims=500, threads=1 )
      # print(bs_cpl$p)

    }
    
    ## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
    par(las=1)
    
    # myboxplot( (s0-s1b)*1e-15 ~ icont, data=df_impacts, ylab="difference s1b - s0 (PgC)", col="grey70", names=cont, ylim=c(-0.01,0.1) )
    # points(1:5, (mean_cont$s0-mean_cont$s1b)*1e-15, pch=16, col="red")
    # abline( h=0, lty=3 )
      
    myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="ratio s1b/s0", col="grey70", names=cont, ylim=c(0.2,3) )
    abline( h=1, lty=3 )
    # points(1:5, mean_cont$ampl_s1b, pch=16, col="red")

    mtext( paste0( letters[icont+1], ")"), font=2, adj = 0, line = 0.5, cex = 1 )
    
  dev.off()


## Plot PDF of x>X
pdf( "fig/extremes.pdf", width=10, height = 7 )
par( mfrow=c(2,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
for (k in 1:nconts) {
# par( mfrow=c(1,3), las=1, mar=c(4,4.5,3,1), mgp=c(3,1,0) )
# for (k in c(1,3,5)) {

  # df <- tibble( prob=(1:n)/sum(1:n), 
  #   s1a=sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, 
  #   s1b=sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, 
  #   s1c=sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9
  #   )
  n0 <- length(IMPACT_s0[[k]])
  n1 <- length(IMPACT_s1c[[k]])
  n <- min(n0,n1)
  plot(   sort( -1 * IMPACT_s1b[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy",ylab="p(x)",xlab="PgC", pch=16, col = rgb(1,0,0,1), axes=FALSE )

  points(   sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )
  points(   sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), log="xy", pch=16, col = rgb(1,0,0,0.4), axes=FALSE )

  # polygon( c( sort( -1 * IMPACT_s1a[[k]][1:n], decreasing=TRUE ) * 1e-9, rev(sort( -1 * IMPACT_s1c[[k]][1:n], decreasing=TRUE ) * 1e-9) ), c( (1:n)/sum(1:n), rev((1:n)/sum(1:n)) ), border = NA, col = rgb(1,0,0,0.3) )
  points( sort( -1 * IMPACT_s0[[k]][1:n],  decreasing=TRUE ) * 1e-9, (1:n)/sum(1:n), pch=16, col = rgb(0,0,0,1) )
  mtext( bquote( alpha["s1b"] == .(format( fit_s1b[[k]]$alpha, digits = 3) ) ), side=1, line=-4, adj=0.1 )
  mtext( bquote( alpha["s0"]  == .(format( fit_s0[[k]]$alpha,  digits = 3) ) ), side=1, line=-2, adj=0.1 )
  mtext( continent[k], line=1, font=2, adj=0 )
  axis(1, mgp=c(3,1,0) )
  axis(2, mgp=c(3.5,1,0))
  box()
}
legend( "left", pch=16, col=c(rgb(1,0,0,1), rgb(0,0,0,1)), legend=c("s1","s0"), bty="n", cex = 1.5 )
dev.off()


## Plot difference s1 - s0, absolute
pdf("fig/extremes_s1-s0_abs_glob.pdf",width=10)
  # par(las=1)
  # with( df_impacts, plot( -1e-15*(s1b-s0), ylab="s1b-s0 (PgC)", xlab="event rank", pch=16, col = rgb(0,0,0,1)) )
  # abline( h=0, lty=3 )
  par(las=1)
  myboxplot( (s0-s1b) ~ icont, data=df_impacts, ylab="difference s1b - s0 (PgC)", col="grey70", names=cont )
  abline( h=0, lty=3 )
  points(1:5, mean_cont$s0-mean_cont$s1b, pch=16, col="red")
dev.off()

## Plot difference s1 - s0, relative
pdf("fig/extremes_s1-s0_rel_glob.pdf",width=10)
  par(las=1)
  with( df_impacts, plot( (s1b/s0), ylab="ratio s1b/s0", xlab="event rank", pch=16, col = rgb(0,0,0,1)) )
  abline( h=1, lty=3 )
dev.off()

## Plot difference s1 - s0, relative (no dependence of amplification factor on size found)
pdf( "fig/extremes_s1-s0_rel_glob.pdf", width=7 )
  par(las=1)
  myboxplot( ampl_s1b ~ icont, data=df_impacts, ylab="ratio s1b/s0", col="grey70", names=cont, ylim=c(0.5,2.5) )
  abline( h=1, lty=3 )
  points(1:5, mean_cont$ampl_s1b, pch=16, col="red")
dev.off()




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
