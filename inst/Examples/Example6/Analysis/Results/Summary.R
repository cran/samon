# -----------------------------------------------------------
options( width=150, max.print=1000 )
library(samon, lib.loc="../../../../samlib")

trt1Results <- readRDS("RDS/treatment1Results.rds")
trt2Results <- readRDS("RDS/treatment2Results.rds")

TM1 <- samonSummaryIM(trt1Results)
TM2 <- samonSummaryIM(trt2Results)

saveRDS(TM1, "RDS/TM1.rds" )
saveRDS(TM2, "RDS/TM2.rds" )

DM <- samonDifferenceSummaryIM(TM1,TM2)
CM <- samonCrossSummaryIM(TM1,TM2)

saveRDS(DM,  "RDS/DM.rds"  )
saveRDS(CM,  "RDS/CM.rds"  )





