## PIM.R IM version
## test samonnGenIMIF, samon_ngenIMIF
## 0. initiate env
## 1. fill in data 
## 2. Popt, Qopt and IF for input data  n times
## 3. Generate NSamples
## 3a Fill in each sample n times and
## 3b. Popt, Qopt and IF 
# -----------------------------------------------------------
options( width=150, max.print=650 )

# Retrieve the SGE_TASK_ID and store it as numeric in SGETID
# -----------------------------------------------------------
 SGETID <- Sys.getenv(c("SGE_TASK_ID"))
 SGENID <- as.numeric(SGETID)

 if ( SGENID <= 100 ) {
     dtype <- 1
     part  <- SGENID
 } else {
     dtype <- 2
     part  <- SGENID - 100
 }

 oname    <- sprintf("RDS/Results_%d_%05d.rds",dtype,part)

# -----------------------------------------------------------
library(samon, lib.loc="../../../samlib")

if ( dtype == 1 ) {
    data("DepWork1")
    data  <- DepWork1
    seed0 <- 2121 + part
    seed1 <- 281
} else if ( dtype == 2 ) {
    data("DepWork2")
    data  <- DepWork2
    seed0 <- 3131 + part
    seed1 <- 427
}

NT           <- ncol(data)
inmodel      <- matrix(1,NT,6)
inmodel[1,]  <- 0
inmodel[NT,] <- 0
inmodel[NT-1,4:6] <- 0

alphaList <- -10:10

Results <- samonIM(
    mat             =      data,   # imput matrix
    Npart           =        10,   # number of partitions
    
    InitialSigmaH   =      0.31,   # initial value
    HighSigmaH      =      10.0,   # high value for H
    
    InitialSigmaF   =      0.31,   # initial value
    HighSigmaF      =      10.0,   # high value for F
    
    lb              =         0,   # parameters for
    ub              =         1,   # cumulative 
    zeta1           =       1.0,   # beta distribution
    zeta2           =       1.0,

    NSamples        =        20,   # bootstraps
    NIMimpute       =         5,
    
    seed0           =     seed0,
    seed1           =     seed1,
    inmodel         =   inmodel,   # input model

    retIFiM         =     FALSE,
    retIFiS         =     FALSE,
    
    alphaList       = alphaList
                     )
# -----------------------------------------------------------
print(Results)
saveRDS(Results, file=oname)
