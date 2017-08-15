# -------------------------------------------------------
# Put together results from various runs of samon.
# -------------------------------------------------------
options( width=150, max.print=250 )

library(samon, lib.loc="../../../../samlib")

## there are 100 files each with 20 by 6 imputations
filenames1 <- sprintf("../RDS/Results_1_%05d.rds",1:100)
filenames2 <- sprintf("../RDS/Results_2_%05d.rds",1:100)

trt1 <- samonCombineIM( filenames1 )
trt2 <- samonCombineIM( filenames2 )

saveRDS( trt1, file = "RDS/treatment1Results.rds")
saveRDS( trt2, file = "RDS/treatment2Results.rds")
