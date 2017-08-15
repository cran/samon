# --------------------------------------------------------------------------------
# Takes a list of rds filenames and combines the objects.
# The first file contains the base object and subsequent objects are are checked
# for consistency against this. 
# Bootstrap sample numbers are incremented so that they are unique in the output
# object.
#
# Returns a list.
#
# The following items are checked for consistency across files:
#  n0              -- number of rows in input matrix mat
#  nt              -- number of columns in input matrix mat
#  mat             -- input matrix itself
#  lb              -- lower bound for mapping Y to the [0,1] interval
#  ub              -- upper bound for mapping Y to the [0,1] interval
#  zeta1           -- first parameter to cumulative beta function
#  zeta2           -- second parameter to cumulative beta function
#  Npart           -- Number of partitions
#
#  Should any of these differ between the base file and one of the other files
#  data will still be appended if ForceAppend = TRUE
#
#  The following items are also cheked for consistency but appending is done
#  even if they don't match. A warning is issued.
#  MaxIter         -- maximum number of iterations in Newton's method
#  FAconvg         -- convergence criterion:  absolute change in function
#  FRconvg         -- convergence criterion:  relative change in function
#  SAconvg         -- convergence criterion:  absolute change in slope
#  InitialSigmaH   -- initial value of sigmaH
#  HighSigmaH      -- highest value of sigmaH to consider
#  InitialSigmaF   -- initial value of sigmaF
#  HighSigmaF      -- highest value of sigmaF to consider
#  Nalpha          -- number of alphas
#  alphaList       -- vector of length Nalpha of alphas
#
#  The following items are taken from the base file.  These items may exist
#  in other files, but are ignored:
#  HM              -- optimal sigmaH for mat
#  FM              -- optimal sigmaF for mat
#  IFM             -- IF estimates for mat
#  retIFiM         -- are there IFiM estimates
#
#  if as.numeric(retIFiM) is 1 then the following object is taken from the
#  base file.
#  IFiM            -- individual IF estimates for mat
#
#  all files are checked to see if there are bootstrap sample results.
#  This is achieved by examining the NSamples item to see if it is greater
#  than 0.  Should it be so then the following items are appended to those
#  found in the base file:
#  NSamples        -- number of samples requested 
#  seed0           -- input seed used
#  HS              -- optimal sigmaH for bootstrap samples
#  FS              -- optimal sigmaF for bootstrap samples
#  IFS             -- IF estimates for bootstrap samples
#
#  if as.numeric(retIFiS) is 1 then return IFis
#  retIFiS         -- indicator if individual IF values are returned for
#                     each bootstrap sample
#  IFiS            -- individual IF values for each bootstrap sample
# --------------------------------------------------------------------------------
samonCombineIM <- function( filenames, replaceSampleNo=TRUE, ForceAppend = FALSE, objlist = NULL ) {

  inlist <- (!is.null(objlist))

  if ( inlist == 1 ) {
      n <- length(objlist)
  } else {
      n <- length(filenames)
  }
  
  for ( i in 1:n ) {
     
  if ( inlist == 1 ) {
      Res <- objlist[[i]]
  } else {
      Res <- readRDS(filenames[i])
  }

     if ( i == 1 ) {

       Retobj      <- Res
       objnames0   <- names(Res)

       check00     <- c( Res$n0, Res$nt, Res$lb, Res$ub, Res$zeta1, Res$zeta2, Res$Npart )
       checkmat0   <- Res$mat
       check01     <- c( Res$MaxIter,Res$FAconvg, Res$FRconvg, Res$SAconvg, Res$InitialSigmaH,
                         Res$HighSigmaH, Res$InitialSigmaF, Res$HighSigmaF, Res$Nalpha )
       chkalpha0   <-    Res$alphaList

       NSamples0   <- Res$NSamples
       if ( is.null( Res$retIFiM ) ) retIFiM0 <- 0
       else retIFiM0 <- Res$retIFiM
       if ( is.null( Res$retIFiS ) ) retIFiS0 <- 0
       else retIFiS0 <- Res$retIFiS

       RegSjk0Ind <- 0 
       if ( match("RegSjk", objnames0, -1 ) > 0 ) {
           RegSjk     <- Res$RegSjk
           RegSjk0Ind <- 1
       }

       if ( NSamples0 > 0 ) {

         HS      <- Res$HS
         FS      <- Res$FS
         IFS     <- Res$IFS

         if ( retIFiS0 == 1 ) IFiS <- Res$IFiS
         
         oset    <- max( HS[,1] )
         if ( oset == 0 ) oset <- 1
       }
     } else {

       objnames    <- names(Res)
       
       ## check files are the same  
       check10     <- c( Res$n0, Res$nt, Res$lb, Res$ub, Res$zeta1, Res$zeta2, Res$Npart )
       checkmat1   <- Res$mat
       check11     <- c( Res$MaxIter,Res$FAconvg, Res$FRconvg, Res$SAconvg, Res$InitialSigmaH,
                         Res$HighSigmaH, Res$InitialSigmaF, Res$HighSigmaF, Res$Nalpha )
       chkalpha1   <-    Res$alphaList

       checkV0     <- identical(   check10, check00   )
       checkMat    <- identical( checkmat1, checkmat0 )
       checkV1     <- identical(   check11, check01   )
       checkV2     <- identical( chkalpha1, chkalpha0 )

       if ( !checkV0 || !checkMat ) {
           ErrorMessage <- paste("Error: in samonAppend. File ",i," does not match base file")
           warning(ErrorMessage)
           if ( ForceAppend == FALSE ) return(NULL)
       }     
       if ( !checkV1  || !checkV2 ) {
           WarningMessage <- paste("Warning: in samonAppend. File ",i," does not match base file")
           warning(ErrorMessage)
       }
       
       NSamples   <- Res$NSamples
       if ( is.null( Res$retIFiS ) ) retIFiS <- 0
       else retIFiS <- Res$retIFiS

       if ( NSamples0 > 0 && NSamples > 0 ) {

         HSi    <- Res$HS
         FSi    <- Res$FS
         IFSi   <- Res$IFS

         if ( replaceSampleNo == TRUE ) {
             HSi[ ,1]   <- HSi[, 1] + oset
             FSi[ ,1]   <- FSi[, 1] + oset
             IFSi[,1]   <- IFSi[,1] + oset
         }
         HS     <- rbind(HS,HSi)
         FS     <- rbind(FS,FSi)
         IFS    <- rbind(IFS,IFSi)

         if ( retIFiS0 == 1 && retIFiS == 1 ) {
             IFiSi <- Res$IFiS
             if ( replaceSampleNo == TRUE ) IFiSi[,1] <- IFiSi[,1] + oset
             IFiS  <- rbind(IFiS,IFiSi)
         }    
         
       oset        <- max( HSi[,1] )
     }
   }
 }
 if ( NSamples0 > 0 ) { 
  Retobj$HS    <- HS
  Retobj$FS    <- FS
  Retobj$IFS   <- IFS

  if ( retIFiS == 1 ) {
    Retobj$IFiS  <- IFiS
  }
 } 
  return(Retobj)
}
