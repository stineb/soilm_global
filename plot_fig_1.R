load( "data/data_aligned_agg.Rdata" )
source("plot_aligned_all.R")
source("plot_bias_resolved_fLUE.R")

out <- get_bias_resolved_flue( filter( df_dday_8d_agg, mysitename!="US-Var" ) )

magn <- 4
ncols <- 2
nrows <- 1
widths <- c(1.2,1.6)*magn
heights <- rep(magn,nrows)
order <- matrix( c(1,2), nrows, ncols, byrow=FALSE)

pdf( "fig/fig_1.pdf", width=sum(widths), height=sum(heights) )

  panel <- layout(
            order,
            widths=widths,
            heights=heights,
            TRUE
            )
  # layout.show( panel )

  plot_aligned_all( df_dday_agg, df_dday_8d_agg, filn=NA )
  mtext( "a)", font=2, adj = 0, line = 0.5, cex = 1 )

  plot_bias_resolved( out$tmp, out$tmp0, out$tmp1, out$tmp4, out$tmp3, cex=0.8, filn=NA )
  mtext( "b)", font=2, adj = 0, line = 0.5, cex = 1 )

dev.off()