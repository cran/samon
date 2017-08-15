## ------------------------------------------------
## Check the VAS data
## Should have intermittent missing data
## ------------------------------------------------
options( width=150, max.print=1000)

library(samon, lib.loc="../../../samlib")
data("VAS1")
data("VAS2")

# Check data
chk1  <- samonDataCheck( VAS1  )
chk2  <- samonDataCheck( VAS2  )

NT <- ncol(VAS1)
means1   <- matrix(NA, NT, NT )
means2   <- matrix(NA, NT, NT )

for ( l in 1:NT ) {
for ( t in 1:l ) {
  if ( length( chk1$desc[,2]  == l ) > 0 )  means1[l,t]  <- mean( VAS1[ chk1$desc[,2]  == l & !is.na(VAS1[,t]), t ] )  
  if ( length( chk2$desc[,2]  == l ) > 0 )  means2[l,t]  <- mean( VAS2[ chk2$desc[,2]  == l & !is.na(VAS2[,t]), t ] )  
}
}

cols <- c( "#320BD6FF",  "#A18B3DFF", "#EB7D54FF", "#17A645FF", "#E46F90FF", "#9F4CDAFF", "#5016C3FF", "#974654FF", "#6072D2FF")

fn <- function(means,chk,NT,fname,title) {
  pdf(file=fname, height=5.1, width=5.5)
  par(mar=c(4,5,2,0.6))

  plot.new()
  plot.window( xlim = c(1-1,NT-1), ylim = c(25,75) )

  for ( i in 1:NT ) {
      lines(  x = 0:(i-1), y = means[i,1:i], lwd=3,  col = cols[i])
      points( x = 0:(i-1), y = means[i,1:i], pch=20, col = cols[i], cex=1.5)
      
      Nvals <- paste("(",sum(chk$desc[,2] == i),")",sep="")
      if ( i == NT ) {
          text(   x = i - 1 - 0.1, y = means[i,  i], Nvals,col = cols[i], pos=1, cex=1.2)
      }
      else {
          text(   x = i - 1,       y = means[i,  i], Nvals,col = cols[i], pos=4, cex=1.2)
      }
  }
  
  axis(1,cex.axis=1.3)
  axis(2,cex.axis=1.3)

  title( main = title, xlab = "Visit",  ylab = "Mean VAS by Last Observation", cex.lab=1.4, cex.main=1.4, font.main=1)
  box()

  dev.off()
}
##options(cex=1.3)
fn(means1, chk1, NT, "Plots/DataCheck1.pdf", "Treatment 1")
fn(means2, chk2, NT, "Plots/DataCheck2.pdf", "Treatment 2")

