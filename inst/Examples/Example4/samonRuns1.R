# samonRuns1.R
# Produce 50 bootstraps each with jackknife estimates.
# Both the output filename and the seed used are based on the 
# SUN grid engine task id.
# ----------------------------------------------------------------
options( width=150, max.print=1000000 )

# some seeds to make more seeds
Seeds1 <- c(  826847827,  679365679,  779140483, 146694967,  539272129 )

# Retrieve the SGE_TASK_ID and store it as numeric in SGETID
# --------------------------------------

 SGETID   <- Sys.getenv(c("SGE_TASK_ID"))
 SGENID   <- as.numeric(SGETID)

 samp     <- SGENID
 n        <- floor(log(samp,10))+1
 oname    <- paste("RDS/Sample1/results_", paste(rep(0,5-n),sep="",collapse=""), samp, ".rds",sep="")

 sd       <- Seeds1[1] 
 set.seed(sd)
 seedList <- ceiling( 1000000 * runif( 1000 ) )
 seed     <- seedList[SGENID]

# only do the jackknife on the main data once
MJK <- ( SGENID == 1 )

# --------------------------------------
library(samon, lib.loc="../../samlib")
data("samonPANSS1")
PANSS1 <- samonPANSS1

alphaList <- -20:20

Results <- samon(
    mat             =    PANSS1,   # imput matrix
    Npart           =        10,   # number of partitions
    
    InitialSigmaH   =      10.0,   # initial value
    HighSigmaH      =      50.0,   # high value for p
    
    InitialSigmaF   =       8.0,   # initial value
    HighSigmaF      =      50.0,   # high value for q
    
    lb              =        30,   # parameters for
    ub              =       210,   # cumulative 
    zeta1           =       4.0,   # beta distribution
    zeta2           =       7.0,

    NSamples        =        50,   # bootstraps
    seed0           =      seed,

    MJackknife      =       MJK,
    SJackknife      =      TRUE,
    
    alphaList       =  alphaList
                     )
# --------------------------------------
saveRDS(Results, file=oname)
