# Samon Example 2. 
# ----------------------------------------------------
# Produce bias corrected estimates: PANSS group 1 and 2.
# ----------------------------------------------------
options( width=150, max.print=1000000 )

library(samon, lib.loc="../../samlib")

# Treatment 1 and 2 data.
data(samonPANSS1)
data(samonPANSS2)

alphaList <- -10:10

Results1 <- samon(mat = samonPANSS1, Npart = 10,
    
    InitialSigmaH   = 10.0,    # initial value
    HighSigmaH      = 50.0,    # high value for H
    
    InitialSigmaF   = 8.0,     # initial value
    HighSigmaF      = 50.0,    # high value for F
    
    lb              = 30,      # parameters to
    ub              = 210,     # cumulative 
    zeta1           = 4.0,     # beta distribution
    zeta2           = 7.0,

    alphaList       = alphaList                  
                     )

Results2 <- samon(mat = samonPANSS2, Npart = 10,
    
    InitialSigmaH   = 10.0,    # initial value
    HighSigmaH      = 50.0,    # high value for H
    
    InitialSigmaF   = 8.0,     # initial value
    HighSigmaF      = 50.0,    # high value for F
    
    lb              = 30,      # parameters to
    ub              = 210,     # cumulative 
    zeta1           = 4.0,     # beta distribution
    zeta2           = 7.0,

    alphaList       = alphaList                  
                     )

print(Results1$IFM)
print(Results2$IFM)
# ----------------------------------------------------

# Select the IF results from results
Trt1IFResults <- Results1$IFM
Trt2IFResults <- Results2$IFM

pdf(file="Plots/Example2.pdf", height=5.5, width=5.7)
par(mar=c(4,5,0.6,0.6))

plot.new()
plot.window( xlim = c(-10, 10), ylim = c( 60, 100 ))

lines( x = Trt1IFResults[,3], y = Trt1IFResults[,6], lwd=3, col = "#CC8866FF")
lines( x = Trt2IFResults[,3], y = Trt2IFResults[,6], lwd=3, col = "#6688CCFF")

axis(1,cex.axis=1.3)
axis(2,cex.axis=1.3)

legend( -9.2, 99, legend = c("treatment 1", "treatment 2"), xjust=0, yjust=1, lty = c("solid","solid"), lwd=c(5,5), col=c("#CC8866FF","#6688CCFF"))

title( xlab = expression(alpha),  ylab = "PANSS Estimate at Visit 5",cex.lab=1.4)
box()

dev.off()
