# -----------------------------------------------------------
options( width=150, max.print=1000 )
library(samon, lib.loc="../../../../samlib")

trt1Results <- readRDS("RDS/treatment1Results.rds")

HM  <- trt1Results$HM
FM  <- trt1Results$FM
IFM <- trt1Results$IFM

HS  <- trt1Results$HS
FS  <- trt1Results$FS
IFS <- trt1Results$IFS

##trt1Results$HM  <- HM [ HM [,2] < 6, ]
##trt1Results$FM  <- FM [ FM [,2] < 6, ]
##trt1Results$IFM <- IFM[ IFM[,2] < 6, ]

##trt1Results$HS  <- HS [ HS [,2] < 6, ]
##trt1Results$FS  <- FS [ FS [,2] < 6, ]
##trt1Results$IFS <- IFS[ IFS[,2] < 6, ]

##trt1Results$NFills <- 5

trt2Results <- readRDS("RDS/treatment2Results.rds")

HM  <- trt2Results$HM
FM  <- trt2Results$FM
IFM <- trt2Results$IFM

HS  <- trt2Results$HS
FS  <- trt2Results$FS
IFS <- trt2Results$IFS

##trt2Results$HM  <- HM [ HM [,2] < 6, ]
##trt2Results$FM  <- FM [ FM [,2] < 6, ]
##trt2Results$IFM <- IFM[ IFM[,2] < 6, ]

##trt2Results$HS  <- HS [ HS [,2] < 6, ]
##trt2Results$FS  <- FS [ FS [,2] < 6, ]
##trt2Results$IFS <- IFS[ IFS[,2] < 6, ]

##trt2Results$NFills <- 5

TM1 <- samonSummaryIM(trt1Results)
TM2 <- samonSummaryIM(trt2Results)

saveRDS(TM1, "RDS/TM1.rds" )
saveRDS(TM2, "RDS/TM2.rds" )

DM <- samonDifferenceSummaryIM(TM1,TM2)
CM <- samonCrossSummaryIM(TM1,TM2)

saveRDS(DM,  "RDS/DM.rds"  )
saveRDS(CM,  "RDS/CM.rds"  )





