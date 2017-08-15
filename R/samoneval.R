samoneval <- function(mat,
         Npart = 10,
         sigmaList = c(1.0),             
         type = "both"
                  )
{
    
  mat    <- as.matrix(mat)
  nt     <- ncol(mat)
  n0     <- nrow(mat)

##  anynaMat <- any( is.na(mat) )  
##  if ( anynaMat ) {
##    SamonMinOffset    <- min(mat,na.rm=TRUE)
##    mat               <- mat - SamonMinOffset + 1
##    mat[ is.na(mat) ] <- -1
##  }

  nsigma <- length(sigmaList)

  noutrows  <- nsigma;

  if ( match( type, c("H","h"), -1 ) > 0 ) {
      Dotype <- 1
  } else if ( match( type, c("F","f"), -1 ) > 0 ) {
      Dotype <- 2
  } else {
      Dotype <- 3
  }
  
  outmat    <- matrix(0,nrow=noutrows,ncol=2)

  if ( Dotype == 1 ) {

    outmat <- matrix(0,nrow=noutrows,ncol=2)
    HFtp   <- 1
      
    fout <-  .C("samon_eval",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          outmat            = as.double  ( as.vector(    outmat ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
          NAOK=TRUE)

    samonOut <- matrix(fout$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    colnames(samonOut) <- c("sigmaH","lossH")
    
  } else if ( Dotype == 2 ) {

    outmat <- matrix(0,nrow=noutrows,ncol=2)
    HFtp   <- 2
      
    fout <-  .C("samon_eval",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          outmat            = as.double  ( as.vector(    outmat ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),      
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
          NAOK=TRUE)

    samonOut <- matrix(fout$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    colnames(samonOut) <- c("sigmaF","lossF")
    
   } else  if ( Dotype == 3 ) {

    outmatH <- matrix(0,nrow=noutrows,ncol=2)
    HFtp    <- 1
      
    foutH <-  .C("samon_eval",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          outmat            = as.double  ( as.vector(   outmatH ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),      
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
          NAOK=TRUE)

    outmatF <- matrix(0,nrow=noutrows,ncol=2)
    HFtp    <- 2
      
    foutF <-  .C("samon_eval",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          outmat            = as.double  ( as.vector(   outmatF ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),      
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
          NAOK=TRUE)

    outmatH  <- matrix(foutH$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    outmatF  <- matrix(foutF$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    samonOut <- cbind(outmatH,outmatF)
    colnames(samonOut) <- c("sigmaH","lossH","sigmaF","lossF")
   }
  return(samonOut)
}
