options( width=150, max.print=1000000 )

library(samon, lib.loc="../../samlib")
data("samonPANSS1")
data("samonPANSS2")

# original data
Y1  <- as.matrix(samonPANSS1)
Y2  <- as.matrix(samonPANSS2)

# samon results
samonR1  <- readRDS("RDS/TM1.rds")
samonR2  <- readRDS("RDS/TM2.rds")
samonRD  <- readRDS("RDS/DM.rds")

alpha0Results <- alpha0TableResults( Y1, samonR1, Y2, samonR2, samonRD )
alpha0Tables(alpha0Results)

    
