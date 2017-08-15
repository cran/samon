options( width=150, max.print=1000000 )

library(samon, lib.loc="../../samlib")
data("samonPANSS1")
data("samonPANSS2")

# original data
Y1  <- as.matrix(samonPANSS1)
Y2  <- as.matrix(samonPANSS2)

TM1 <- readRDS( "RDS/TM1.rds" )
TM2 <- readRDS( "RDS/TM2.rds" )

TM1 <- TM1$TM
TM2 <- TM2$TM

Results1 <- TM1[,c("alpha","IFEst")]
Results2 <- TM2[,c("alpha","IFEst")]

difference1 <- samonECompleterStatus( Y1[,ncol(Y1)], Results1[,2] )
difference2 <- samonECompleterStatus( Y2[,ncol(Y2)], Results2[,2] )

# ----------------------------------------------------------------------------------------
pdf(file="Plots/EYdiffCompleterStatus.pdf", height=5.5, width=6)
par(mar=c(4,5,1,1))

plot.new()
plot.window( xlim = range(Results1[,1]), ylim = range(c(difference1, difference2)) )

lines( x = Results1[,1], y = difference1, lwd=3, col = "#CC8866FF")
lines( x = Results2[,1], y = difference2, lwd=3, col = "#6688CCFF")

legend( -19, max( c( difference1, difference2 ) ) - 0.1,
       legend = c("Placebo", "Active"),
       xjust=0, yjust=1, lty = c("solid","solid"), lwd=c(5,5,5), col=c("#CC8866FF","#6688CCFF"))

axis(1)
axis(2)

title( xlab = expression(alpha),
       ylab = "Difference in Means\n(Non-completers minus Completers)")
box()
dev.off()
