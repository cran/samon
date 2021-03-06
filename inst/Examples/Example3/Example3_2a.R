# Example3_2a.R
# Produce bias corrected influence function estimates
# and 500 bootstrap estimates for treatment 2.
# ----------------------------------------------------
options( width=150, max.print=10000 )

library(samon, lib.loc="../../samlib")

alphaList <- -10:10

# Treatment 2 data.
data(samonPANSS2)

Results2a <- samon(
    mat             =  samonPANSS2,   # imput matrix
    Npart           =  10,            # number of partitions
    
    InitialSigmaH   =  10.0,          # initial value
    HighSigmaH      =  50.0,          # high value for H
    
    InitialSigmaF   =  8.0,           # initial value
    HighSigmaF      =  50.0,          # high value for F
    
    lb              =  30,            # parameters to
    ub              =  210,           # cumulative 
    zeta1           =   4.0,          # beta distribution
    zeta2           =   7.0,

    NSamples        =     500,        # bootstraps
    seed0           =  281881,

    MJackknife      =     TRUE,       # jackknife main data
    SJackknife      =     TRUE,       # jackknife bootstraps
    retIFiM         =    FALSE,
    retIFiS         =    FALSE,
    
    alphaList       =  alphaList
    
                     )
# --------------------------------------

# save results for later use
saveRDS(Results2a, "RDS/Results2a.rds")
