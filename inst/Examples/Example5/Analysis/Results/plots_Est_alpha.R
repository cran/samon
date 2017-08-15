# --------------------------------------
options( width=150, max.print=250 )

TM1 <- readRDS("RDS/TM1.rds"  )
TM2 <- readRDS("RDS/TM2.rds"  )
DM  <- readRDS("RDS/DM.rds"   )

CI1 <- TM1$CI
TM1 <- TM1$TM
CI2 <- TM2$CI
TM2 <- TM2$TM
CID <- DM$CI
TMD <- DM$TM

# alpha estimates and confidence intervals
Results1 <- cbind( TM1[,c("alpha","IFEst")     ],  CI1[,c("lb7","ub7")] )
Results2 <- cbind( TM2[,c("alpha","IFEst")     ],  CI2[,c("lb7","ub7")] )
ResultsD <- cbind( TMD[,c("alpha","Difference")],  CID[,c("lb7","ub7")] )

# Takes a matrix of results, Res, with Res[,1] to be plotted on the xaxis,
# Res[,2] (estimate) to be plotted on the yaxis.  Res[,6] and Res[,7]
# contain lower and upper bounds to be plotted as a band on the plot.
samonPlot <- function( Res, file, height, width, ylab, xlim, ylim, legpos, maintext) {
    
  pdf(file=file, height=height, width=width)
  par(mar=c(4,5,0.6,0.6))

  plot.new()
  plot.window( xlim = xlim, ylim = ylim )

  lines( x = Res[,1], y = Res[,3], lwd=3, lty = c("solid"), col = "#77AAFFFF")
  lines( x = Res[,1], y = Res[,4], lwd=3, lty = c("solid"), col = "#77AAFFFF")

  lines( x = Res[,1], y = Res[,2], lwd=3, col = "#CC8866FF")

  axis(1,cex.axis=1.4)
  axis(2,cex.axis=1.4)

  title( main = maintext, xlab = expression(alpha),  ylab = ylab, cex.lab=1.5)
  box()

  dev.off()
  invisible(return())
}

samonPlot(Results1, "Plots/Results1.pdf",  5.5, 6, "Estimate",    c(-10,10), c(  20.0,  60.0 ), c( 3.6,  0.85 ), maintext = NULL )
samonPlot(Results2, "Plots/Results2.pdf",  5.5, 6, "Estimate",    c(-10,10), c(  20.0,  60.0 ), c( 3.6,  0.85 ), maintext = NULL )
samonPlot(ResultsD, "Plots/ResultsD.pdf",  5.5, 6, "Difference",  c(-10,10), c( -20.0,  10.0), c( 3.6, -0.05), maintext = NULL )

