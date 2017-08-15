# --------------------------------------------------------------------------------
# Plot the P and Q smoothing functions for a range of smoothing parameters.
# --------------------------------------------------------------------------------
options( width=150, max.print=500 )


# gather the smoothing results.
fn1 <- function( i ) {
    ifile <- sprintf("RDS/HFResults_1.%05d.rds", i )
    trtR  <- readRDS(ifile)
    return( trtR )
}

fn2 <- function( i ) {
    ifile <- sprintf("RDS/HFResults_2.%05d.rds", i )
    trtR  <- readRDS(ifile)
    return( trtR )
}

trt1R <- lapply(1:5,fn1)
trt2R <- lapply(1:5,fn2)

trt1R <- do.call(rbind,trt1R)
trt2R <- do.call(rbind,trt2R)

print( trt1R )
print( trt2R )

colnames(trt1R) <- c("impute", "dtype", "sigmaH", "lossH", "sigmaF", "lossF")
colnames(trt2R) <- c("impute", "dtype", "sigmaH", "lossH", "sigmaF", "lossF")

ag.o <- aggregate( lossH ~ impute + dtype, trt1R, min )
print(ag.o)
Hmin1 <- merge(trt1R,ag.o,sort=FALSE)
Hmin1 <- Hmin1[,c("impute","dtype","sigmaH","lossH")]
print(Hmin1)

ag.o <- aggregate( lossF ~ impute + dtype, trt1R, min )
print(ag.o)
Fmin1 <- merge(trt1R,ag.o,sort=FALSE)
Fmin1 <- Fmin1[,c("impute","dtype","sigmaF","lossF")]
print(Fmin1)

ag.o <- aggregate( lossH ~ impute + dtype, trt2R, min )
print(ag.o)
Hmin2 <- merge(trt2R,ag.o,sort=FALSE)
Hmin2 <- Hmin2[,c("impute","dtype","sigmaH","lossH")]
print(Hmin2)

ag.o <- aggregate( lossF ~ impute + dtype, trt2R, min )
print(ag.o)
Fmin2 <- merge(trt2R,ag.o,sort=FALSE)
Fmin2 <- Fmin2[,c("impute","dtype","sigmaF","lossF")]
print(Fmin2)
                              


cols <- c("#ffb3baFF","#ffdfbaFF","#ffffbaFF","#baffc9FF","#bae1ffFF","#edc951FF","#eb6841FF","#cc2a36FF","#4f372dFF","#00a0b0FF",
          "#757676FF","#800909FF","#e0cda7FF","#2a334fFF","#ac8f57FF","#00b159FF","#005b96FF","#6497b1FF","#a32020FF","#602320FF",
          "#666547FF","#ffe28aFF","#fffeb3FF","#fb2e01FF","#6fcb9fFF","#e8d174FF","#e39e54FF","#d64d4dFF","#4d7358FF","#9ed670FF",
          "#a67c00FF","#ffbf00FF")

cols <- c("#415e3cFF", "#8d0808FF", "#0f50a0FF", "#007fccFF", "#3c415aFF", "#fff68fFF", "#8a2be2FF" )
cols <- c("#585858FF", "#118C4EFF", "#C1E1A6FF", "#FF9009FF", "#6DBDD6FF", "#B71427FF", "#FFE658FF" )

HFPlot <- function( Res, file, height, width, ylab, xlim, ylim, legpos, nimputes, HFmin, rc) {
    
  pdf(file=file, height=height, width=width)
  par(mar=c(4,5,0.6,0.6))

  plot.new()
  plot.window( xlim = xlim, ylim = ylim )

  axis(1,cex.axis=1.2)
  axis(2,cex.axis=1.2)

  for ( i in 1:nimputes ) {
    XY <- Res[ Res[,1] == i, ]  
    lines( x = XY[,2], y = XY[,3], lwd=1, col = cols[(i%%rc)+1], lty=c("solid"))
    rug( HFmin[i,3],col=cols[(i%%rc)+1])
  }
  title( xlab = expression(sigma),  ylab = ylab, cex.lab=1.2)
  box()

  dev.off()
  invisible(return())
}

HFPlot(trt1R[,c(1,3,4)], "Plots/H1.pdf", 4.7, 5.0, "Loss Function (H)",  c(0,5), c( 1.8845, 1.89), c(25, 8), 5, Hmin1, rc = 17 )
HFPlot(trt1R[,c(1,5,6)], "Plots/F1.pdf", 4.7, 5.0, "Loss Function (F)",  c(0,5), c( 2.49, 2.7), c(25,11), 5, Fmin1, rc = 17 )

HFPlot(trt2R[,c(1,3,4)], "Plots/H2.pdf", 4.7, 5.0, "Loss Function (H)",  c(0,5), c( 1.973, 1.976), c(25, 8), 5, Hmin2, rc = 17 )
HFPlot(trt2R[,c(1,5,6)], "Plots/F2.pdf", 4.7, 5.0, "Loss Function (F)",  c(0,5), c( 2.16,  2.32), c(25,11), 5, Fmin2, rc = 17 )

