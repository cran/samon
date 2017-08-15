# ----------------------------------------------------------------
# Produces estimates from samon output matrices for treatment1
# treatment2 and their difference treatment2 - treatment1.
# ----------------------------------------------------------------

options( width=200, max.print=700 )

library(samon,lib.loc="../../samlib")

# treatment 1, 2 and difference
# ---------------------------------------------------
trt1Results <- readRDS("RDS/treatment1Results.rds")
trt2Results <- readRDS("RDS/treatment2Results.rds")

TM1 <- samonSummary(trt1Results)
TM2 <- samonSummary(trt2Results)

DM  <- samonDifferenceSummary(TM1,TM2)
CM  <- samonCrossSummary(TM1,TM2)

saveRDS(TM1, "RDS/TM1.rds" )
saveRDS(TM2, "RDS/TM2.rds" )

saveRDS(DM,  "RDS/DM.rds"  )
saveRDS(CM,  "RDS/CM.rds"  )

