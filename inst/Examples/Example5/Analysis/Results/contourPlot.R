options( width=200, max.print=10000 )

# --------------------------------------

XRes <- readRDS("RDS/CM.rds")
XRes <- as.matrix(XRes$CI)

print(XRes)

pdf(file="Plots/contour.pdf", height=5, width=6.8)
par(mar=c(4,5,0.5,0.5),cex.axis=1.2,cex.lab=1.4)

filled.contour( x = -10:10,
                y = -10:10,
                z = matrix(XRes[,"Difference"],21,21, byrow=TRUE),
                xlab    = expression(paste(alpha, " (Placebo)")),
                ylab    = expression(paste(alpha, " (Active)")),
                nlevels = 8,
                color.palette     = colorRampPalette(c( "#993404","#D95F0E","#FE9929","#FFD9BE","#FFFFD4"), space="rgb"),
                plot.axis = ( points( XRes[ sign(XRes[,"lb7"]) == sign(XRes[,"ub7"]), c(1,2)], pch=15, cex=0.6, col = c("#44447799")))
              ) 
dev.off()


