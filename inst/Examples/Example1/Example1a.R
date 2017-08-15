# Samon Example 1a. 
# --------------------------------------------------------------
# Examining the loss functions as functions of their
# smoothing parameters.
# --------------------------------------------------------------
options( width=150, max.print=1000000 )

library(samon, lib.loc="../../samlib")

HFPlot <- function( Res, file, height, width, ylab, xlim, ylim, legpos) {
    
  pdf(file=file, height=height, width=width)
  par(mar=c(4,5,0.6,0.6))

  plot.new()
  plot.window( xlim = xlim, ylim = ylim )

  axis(1,cex.axis=1.2)
  axis(2,cex.axis=1.2)

  lines( x = Res[,1], y = Res[,2], lwd=3, col = "#EE8855FF")
  lines( x = Res[,3], y = Res[,4], lwd=3, col = "#888888FF")

  legend( legpos[1], legpos[2],
          legend = c("Treatment 1", "Treatment 2"),
          xjust=0, yjust=1, lty = c("solid","solid"),
          lwd=c(5,5), col=c("#EE8855FF","#888888FF"))
  
  title( xlab = expression(sigma),  ylab = ylab, cex.lab=1.2)
  box()

  dev.off()
  invisible(return())
}

# --------------------------------------

data("samonPANSS1")
HF1 <- samoneval( mat = samonPANSS1, Npart = 10,
       sigmaList       = seq(0.5,50.0,by=0.1),          
       type            = "both" )

data("samonPANSS2")
HF2 <- samoneval( mat = samonPANSS2, Npart = 10,
       sigmaList       = seq(0.5,50.0,by=0.1),          
       type            = "both" )

# --------------------------------------

ResultsH <- cbind(HF1[,c(1,2)],HF2[,c(1,2)])
ResultsF <- cbind(HF1[,c(3,4)],HF2[,c(3,4)])

HFPlot(ResultsH, "Plots/H.pdf", 4.7, 5.0, "Loss Function (H)",  c(0,50), c( 3, 8), c(25,8) )
HFPlot(ResultsF, "Plots/F.pdf", 4.7, 5.0, "Loss Function (F)",  c(0,50), c( 3, 8), c(25,8) )

