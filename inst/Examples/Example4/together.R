# --------------------------------------------------------------------------------
# Put together results from various runs of samon.
#
# Makes  treatment1Results.rds  with treatment 1 results
#        treatment2Results.rds  with treatment 2 results
# --------------------------------------------------------------------------------

options( width=150, max.print=250 )

library(samon, lib.loc="../../samlib")

# there are 100 files each with 50 bootstraps per treatment arm.
filenames1 <- sprintf("RDS/Sample1/results_%05d.rds", 1:100 )
filenames2 <- sprintf("RDS/Sample2/results_%05d.rds", 1:100 )

trt1 <- samonCombine( filenames1 )
trt2 <- samonCombine( filenames2 )

saveRDS( trt1, file = "RDS/treatment1Results.rds")
saveRDS( trt2, file = "RDS/treatment2Results.rds")
