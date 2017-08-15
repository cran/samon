options( width=150, max.print=1000 )

library(samon, lib.loc="../../../samlib")
data("VAS1")
data("VAS2")

Y1 <- VAS1
Y2 <- VAS2

Tmat1 <- samonTabmat1( Y1 )
Tmat2 <- samonTabmat1( Y2 )

samonTable1( Tmat1, trtlab = "Treatment 1" )
samonTable1( Tmat2, trtlab = "Treatment 2" ) 

