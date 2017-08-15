## -----------------------------------------------------
## Produce a summary table for a matrix Y of data.
##
## Usage: samonTabmat1( Y )
##
## Y is an N0 by NT matrix with columns representing
## timepoints and rows representing subjects.
##
## Output matrix is NT by 10 where NT is the number
## of timepoints in a study and NT = ncol(Y)
## 
## The 10 columns have 
##  1. t
##  2. Number On Study
##  3. Number Observed
##  4. last seen
##  5. Proportion last seen (of on-study) 
##  6. Proportion last seen (of observed)
##  7. inter miss
##  8. proportion inter missing (of onstudy)
##  9. mean observed
## 10. std observed
## -----------------------------------------------------
samonTabmat1 <- function( Y ) {
    
  NT <- ncol(Y)
  N0 <- nrow(Y)
  
  L <- rep(0,N0)
  for ( i in 1:NT ) {
    L[ !is.na(Y[,i]) ] <- i
  }  

  MatTab <- matrix( 0, NT, 10 )
  for ( i in 1:NT ) {
    # t  
    MatTab[i,1] <- i   

    # on study
    MatTab[i,2] <- sum( L >= i )

    # values on study at time i
    Vi <- Y[ L >= i, i ]

    # N observed
    MatTab[i, 3] <- sum( !is.na(Vi) )
    
    # last seen
    MatTab[i, 4] <- sum( L == i )
    MatTab[i, 5] <- MatTab[i,4] / MatTab[i,2]
    MatTab[i, 6] <- MatTab[i,4] / MatTab[i,3]

    # intermittent missing
    MatTab[i, 7] <- sum( is.na(Vi) )
    MatTab[i, 8] <- MatTab[i,7] / MatTab[i,2]

    # mean and sd
    MatTab[i, 9] <- mean( Vi[!is.na(Vi)] )
    MatTab[i,10] <- sd  ( Vi[!is.na(Vi)] )
  }
  return(MatTab)
}
