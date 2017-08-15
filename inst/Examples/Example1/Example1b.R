# Samon Example 1b. 
# --------------------------------------------------------------
# Finding optimal smoothing parameters.
# --------------------------------------------------------------
options( width=150, max.print=1000000 )

library(samon, lib.loc="../../samlib")

# Treatment 1
data("samonPANSS1")
print(samonPANSS1)

samonResults <- samon( mat = samonPANSS1, Npart = 10,
    
    InitialSigmaH   = 10.0,    # initial value
    HighSigmaH      = 50.0,    # high value
    
    InitialSigmaF   = 8.0,     # initial value
    HighSigmaF      = 50.0,    # high value

    SAconvg         = 1E-6,    # stopping criteria
    FAconvg         = 1E-6,
    FRconvg         = 1E-6
                      
                     )
# ----------------------------------------------------

print(samonResults$HM)
print(samonResults$FM)

