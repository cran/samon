options( width=150, max.print=1000 )

library(samon, lib.loc="../../../samlib")

data("DepWork1")
data("DepWork2")

Y1 <- DepWork1
Y2 <- DepWork2


Tmat1 <- samonTabmat1( Y1 )
Tmat2 <- samonTabmat1( Y2 )

samonTable1( Tmat1, trtlab = "Treatment 1" )
samonTable1( Tmat2, trtlab = "Treatment 2" ) 

