options( width=150, max.print=1000000 )

library(samon, lib.loc="../../../../samlib")

# the data
data("VAS1")
data("VAS2")
Y1  <- VAS1
Y2  <- VAS2

# samon results
samonR1  <- readRDS("RDS/TM1.rds")
samonR2  <- readRDS("RDS/TM2.rds")
samonRD  <- readRDS("RDS/DM.rds")

alpha0Results <- alpha0TableResultsIM( Y1, samonR1, Y2, samonR2, samonRD )
alpha0TablesIM(alpha0Results)

    
