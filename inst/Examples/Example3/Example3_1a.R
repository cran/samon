# Example3_1a.R
# Produce bias corrected influence function estimates
# and 500 bootstrap estimates for treatment 1.
# ----------------------------------------------------
options( width=150, max.print=10000 )

library(samon, lib.loc="../../samlib")

alphaList <- -10:10

# Treatment 1 data.
data(samonPANSS1)

Results1a <- samon(
    mat             =  samonPANSS1,   # imput matrix
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
    seed0           =   81881,

    MJackknife      =     TRUE,       # jackknife main data
    SJackknife      =     TRUE,       # jackknife bootstraps
    retIFiM         =    FALSE,
    retIFiS         =    FALSE,
    
    alphaList       =  alphaList
    
                     )
# --------------------------------------

# save results for later use
saveRDS(Results1a, "RDS/Results1a.rds")
