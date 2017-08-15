# Example3_plots.R
#
# Results have previously been stored in rds files
# treatment 1:   Results1a.rds
#                Results1b.rds
# treatment 2:   Results2c.rds
#                Results2d.rds
#
# 1. put the results together and plot
# --------------------------------------------------------------------------------
options( width=150, max.print=250 )

library(samon, lib.loc="../../samlib")

# the first two are for treatment 1, the second two for treatment 2.
filenames1 <- c("RDS/Results1a.rds", "RDS/Results1b.rds")
filenames2 <- c("RDS/Results2a.rds", "RDS/Results2b.rds")

trt1Results <- samonCombine( filenames1 )
trt2Results <- samonCombine( filenames2 )

# summarize
Summary1 <- samonSummary( trt1Results )
Summary2 <- samonSummary( trt2Results )
SummaryD <- samonDifferenceSummary( Summary1, Summary2 )

CI1 <- Summary1$CI
TM1 <- Summary1$TM
CI2 <- Summary2$CI
TM2 <- Summary2$TM
CID <- SummaryD$CI
TMD <- SummaryD$TM

# alpha estimates and confidence intervals
Results1 <- cbind( TM1[,c("alpha","IFEst")     ],  CI1[,c("lb8","ub8")] )
Results2 <- cbind( TM2[,c("alpha","IFEst")     ],  CI2[,c("lb8","ub8")] )
ResultsD <- cbind( TMD[,c("alpha","Difference")],  CID[,c("lb8","ub8")] )

# Takes a matrix of results, Res, with Res[,1] to be plotted on the xaxis,
# Res[,2] to be plotted on the yaxis.  Res[,3] and Res[,4] contain lower
# and upper bounds to be plotted as a band on the plot.
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

samonPlot(Results1, "Plots/Trt1Est.pdf",  5.5, 6, "PANSS Estimate (visit 5)",       c(-10,10), c(  60, 110), c( 3.6,  63.0), maintext = NULL )
samonPlot(Results2, "Plots/Trt2Est.pdf",  5.5, 6, "PANSS Estimate (visit 5)",       c(-10,10), c(  60, 110), c( 3.6,  63.0), maintext = NULL )
samonPlot(ResultsD, "Plots/TrtDEst.pdf",  5.5, 6, "Difference in PANSS (visit 5)",  c(-10,10), c( -40,   5), c( 3.6,  -9.0), maintext = NULL )
