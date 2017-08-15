options( width=150, max.print=1000000 )

library(samon, lib.loc="../../../../samlib")

# the data
data("VAS1")
data("VAS2")
Y1  <- VAS1
Y2  <- VAS2

NT   <- ncol(Y1)
Y1NT <- Y1[,NT]
Y2NT <- Y2[,NT]

TM1 <- readRDS( "RDS/TM1.rds" )
TM2 <- readRDS( "RDS/TM2.rds" )

TM1 <- TM1$TM
TM2 <- TM2$TM

Results1 <- TM1[,c("alpha","IFEst")]
Results2 <- TM2[,c("alpha","IFEst")]

difference1 <- samonECompleterStatus( Y1NT, Results1[,2] )
difference2 <- samonECompleterStatus( Y2NT, Results2[,2] )

# ----------------------------------------------------------------------------------------
pdf(file="Plots/EYdiffCompleterStatus.pdf", height=5.5, width=6)
par(mar=c(4,5,1,1))

plot.new()
plot.window( xlim = range(Results1[,1]), ylim = range(c(difference1, difference2)) )

lines( x = Results1[,1], y = difference1, lwd=3, col = "#CC8866FF")
lines( x = Results2[,1], y = difference2, lwd=3, col = "#6688CCFF")

legend( -10 + 0.1, 0.99 * max( c( difference1, difference2 ) ),
       legend = c("Placebo", "Active"),
       xjust=0, yjust=1, lty = c("solid","solid"), lwd=c(5,5,5), col=c("#CC8866FF","#6688CCFF"))

axis(1)
axis(2)

title( xlab = expression(alpha),
       ylab = "Difference in Means\n(Non-completers minus Completers)")
box()
dev.off()
