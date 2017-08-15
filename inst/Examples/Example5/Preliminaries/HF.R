# --------------------------------------------------------------
# Examining the loss functions.
# For a range of smoothing parameters produce plots of loss by
# smoothing parameter.  Multiple imputations are used.
# --------------------------------------------------------------
options( width=150, max.print=1000000 )

# Retrieve the SGE_TASK_ID and store it as numeric in SGETID
# -----------------------------------------------------------
 SGETID <- Sys.getenv(c("SGE_TASK_ID"))
 SGENID <- as.numeric(SGETID)

 oname1 <- sprintf("RDS/HFResults_1.%05d.rds", SGENID )
 oname2 <- sprintf("RDS/HFResults_2.%05d.rds", SGENID )
# -----------------------------------------------------------

library(samon, lib.loc="../../../samlib")

sigmaList <- seq(0.2,40,by=0.1)

# the data
# --------------------------------------
data("VAS1")
data("VAS2")

Y1 <- VAS1
Y2 <- VAS2

NT           <- ncol(Y1)
inmodel      <- matrix(1,NT,6)
inmodel[1,]  <- 0
inmodel[NT,] <- 0
inmodel[NT-1,4:6] <- 0

seeds <- 3121 + 1:100

  seed <- seeds[SGENID]

  HF1 <- samonevalIM( mat = Y1, Npart = 10,
            sigmaList = sigmaList,       
            inmodel = inmodel,
            seed = seed,
            type = "both" )

  out1 <- cbind( SGENID, 1, HF1$OutSig )

  HF2 <- samonevalIM( mat = Y2, Npart = 10,
            sigmaList = sigmaList,       
            inmodel = inmodel,
            seed = seed,
            type = "both" )

  out2 <- cbind( SGENID, 2, HF2$OutSig )

saveRDS(out1,oname1)
saveRDS(out2,oname2)
# --------------------------------------
