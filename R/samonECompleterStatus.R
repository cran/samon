# ----------------------------------------------------------------------------------------
# Expected difference in Y (non-completers - completers).  
samonECompleterStatus <- function( Y, BC ) {
  n0 <- length(Y)

  anynaY <- any( is.na(Y) )  
  if ( anynaY ) {
    ncomplete     <- sum( !is.na(Y) ) 
    nnoncomplete  <- sum(  is.na(Y) )
    EYcomplete    <- mean(Y[ !is.na(Y) ])
  } else {
    ncomplete     <- sum( Y != -1) 
    nnoncomplete  <- sum( Y == -1) 
    EYcomplete    <- mean(Y[ Y != -1 ])
  }
  EYnoncomplete <- ( BC - EYcomplete * ( ncomplete / n0 )) / ( nnoncomplete / n0 )
  difference    <- EYnoncomplete - EYcomplete
  return(difference)
}
