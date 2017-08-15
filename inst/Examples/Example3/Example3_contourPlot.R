# Example3_contourPlot.R
# --------------------------------------
library(samon, lib.loc="../../samlib")

# the first two are for treatment 1, the second two for treatment 2.
filenames1 <- c("RDS/Results1a.rds", "RDS/Results1b.rds")
filenames2 <- c("RDS/Results2a.rds", "RDS/Results2b.rds")

trt1Results <- samonCombine( filenames1 )
trt2Results <- samonCombine( filenames2 )

# summarize
Summary1 <- samonSummary( trt1Results )
Summary2 <- samonSummary( trt2Results )

XRes <- samonCrossSummary( Summary1, Summary2 )
XRes <- as.matrix(XRes$CI)

pdf(file="Plots/Example3_contour.pdf", height=5, width=6.8)
par(mar=c(4,5.0,0.5,0.5),cex.axis=1.2,cex.lab=1.4)

filled.contour(
  x       = -10:10,
  y       = -10:10,
  z       = matrix(XRes[,3],21,21, byrow=TRUE),
  xlab    = expression(paste(alpha, " (Placebo)")),
  ylab    = expression(paste(alpha, " (Active)")),
  nlevels = 8,
  color.palette = colorRampPalette(c( "#993404","#D95F0E","#FE9929",
                                      "#FFD9BE","#FFFFD4"), space="rgb"),
  plot.axis = ( points( XRes[ sign(XRes[,18]) == sign(XRes[,19]), c(1,2)],
                       pch=15, cex=0.6, col = c("#44447799")))
              ) 
dev.off()
