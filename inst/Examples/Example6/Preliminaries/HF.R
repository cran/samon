# --------------------------------------------------------------
# Examining the loss functions.
# For a range of smoothing parameters produce plots of loss by
# smoothing parameter.
# Since the loss function is determined somewhat on the way data
# are filled in, multiple imputations are performed.
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

sigmaList <- seq(0.05,5,by=0.01)

# the data
# --------------------------------------
data("DepWork1")
data("DepWork2")

Y1 <- DepWork1
Y2 <- DepWork2

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
