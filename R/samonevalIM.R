samonevalIM <- function(mat,
         Npart = 10,
         sigmaList = c(1.0),             
         inmodel = inmodel,
         seed = 1,               
         type = "both"
                  )
{
    
  mat <- as.matrix(mat)
  nt  <- ncol(mat)
  n0  <- nrow(mat)

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

  # max model size
  nmodel       <- 6
  FMat         <- matrix(0, n0, nt + 2    )
  LEstsM       <- matrix(0, nt, nmodel + 7)
  ModelsM      <- matrix(0, nt, nmodel    )

  outmat    <- matrix(0,nrow=noutrows,ncol=2)

  if ( Dotype == 1 ) {

    outmat <- matrix(0,nrow=noutrows,ncol=2)
    HFtp   <- 1
      
    fout <-  .C("samon_evalIM",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          nmodel            = as.integer (               nmodel   ),
          inmodel           = as.integer ( as.vector(   inmodel ) ),
          FMat              = as.double  ( as.vector(      FMat ) ),
          LEstsM            = as.double  ( as.vector(    LEstsM ) ),
          ModelsM           = as.integer ( as.vector(   ModelsM ) ),
          outmat            = as.double  ( as.vector(    outmat ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),
          seed              = as.integer (                 seed   ),
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
                NAOK=TRUE)

    OutSig           <- matrix(fout$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    colnames(OutSig) <- c("sigmaH","lossH")
    fout$mat         <- matrix(fout$mat,     nrow= n0,          ncol = nt,       byrow=FALSE  )
    fout$FMat        <- matrix(fout$FMat,    nrow= n0,          ncol = nt + 2,   byrow=TRUE   )
    fout$LEstsM      <- matrix(fout$LEstsM,  nrow= nt,          ncol = nmodel+7, byrow=TRUE   )
    fout$ModelsM     <- matrix(fout$ModelsM, nrow= nt,          ncol = nmodel,   byrow=TRUE   )

    samonOut         <- list( N0 = n0, NT = nt, inmodel = inmodel, OutSig = OutSig, mat = fout$mat, FMat = fout$FMat, LEstsM = fout$LEstsM, ModelsM = fout$ModelsM ) 
  } else if ( Dotype == 2 ) {

    outmat <- matrix(0,nrow=noutrows,ncol=2)
    HFtp   <- 2
      
    fout <-  .C("samon_evalIM",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          nmodel            = as.integer (               nmodel   ),
          inmodel           = as.integer ( as.vector(   inmodel ) ),
          FMat              = as.double  ( as.vector(      FMat ) ),
          LEstsM            = as.double  ( as.vector(    LEstsM ) ),
          ModelsM           = as.integer ( as.vector(   ModelsM ) ),
          outmat            = as.double  ( as.vector(    outmat ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),
          seed              = as.integer (                 seed   ),
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
                NAOK=TRUE)

    OutSig           <- matrix(fout$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    colnames(OutSig) <- c("sigmaF","lossF")
    fout$mat         <- matrix(fout$mat,     nrow= n0,          ncol = nt,       byrow=FALSE  )
    fout$FMat        <- matrix(fout$FMat,    nrow= n0,          ncol = nt + 2,   byrow=TRUE   )
    fout$LEstsM      <- matrix(fout$LEstsM,  nrow= nt,          ncol = nmodel+7, byrow=TRUE   )
    fout$ModelsM     <- matrix(fout$ModelsM, nrow= nt,          ncol = nmodel,   byrow=TRUE   )

    samonOut         <- list( N0 = n0, NT = nt, inmodel = inmodel, OutSig = OutSig, mat = fout$mat, FMat = fout$FMat, LEstsM = fout$LEstsM, ModelsM = fout$ModelsM ) 
   } else  if ( Dotype == 3 ) {

    outmatH <- matrix(0,nrow=noutrows,ncol=2)
    HFtp    <- 1
      
    foutH <-  .C("samon_evalIM",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          nmodel            = as.integer (               nmodel   ),
          inmodel           = as.integer ( as.vector(   inmodel ) ),
          FMat              = as.double  ( as.vector(      FMat ) ),
          LEstsM            = as.double  ( as.vector(    LEstsM ) ),
          ModelsM           = as.integer ( as.vector(   ModelsM ) ),
          outmat            = as.double  ( as.vector(   outmatH ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),
          seed              = as.integer (                 seed   ),
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
                NAOK=TRUE)

    outmatF <- matrix(0,nrow=noutrows,ncol=2)
    HFtp    <- 2
      
    foutF <-  .C("samon_evalIM",
          n0                = as.integer (                   n0   ), 
          nt                = as.integer (                   nt   ),  
          mat               = as.double  ( as.vector(       mat ) ),
          nmodel            = as.integer (               nmodel   ),
          inmodel           = as.integer ( as.vector(   inmodel ) ),
          FMat              = as.double  ( as.vector(      FMat ) ),
          LEstsM            = as.double  ( as.vector(    LEstsM ) ),
          ModelsM           = as.integer ( as.vector(   ModelsM ) ),
          outmat            = as.double  ( as.vector(   outmatF ) ),
          nsigma            = as.integer (               nsigma   ),
          sigmaList         = as.double  ( as.vector( sigmaList ) ),
          seed              = as.integer (                 seed   ),
          Npart             = as.integer (                Npart   ),
          HFtp              = as.integer (                 HFtp   ),
                NAOK=TRUE)
                 

    outmatH  <- matrix(foutH$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    outmatF  <- matrix(foutF$outmat, nrow=noutrows, ncol=2, byrow=TRUE )
    OutSig   <- cbind(outmatH,outmatF)
    colnames(OutSig) <- c("sigmaH","lossH","sigmaF","lossF")

    foutF$mat         <- matrix(foutF$mat,     nrow= n0,          ncol = nt,       byrow=FALSE  )
    foutF$FMat        <- matrix(foutF$FMat,    nrow= n0,          ncol = nt + 2,   byrow=TRUE   )
    foutF$LEstsM      <- matrix(foutF$LEstsM,  nrow= nt,          ncol = nmodel+7, byrow=TRUE   )
    foutF$ModelsM     <- matrix(foutF$ModelsM, nrow= nt,          ncol = nmodel,   byrow=TRUE   )

    samonOut         <- list( N0 = n0, NT = nt, inmodel = inmodel, OutSig = OutSig, mat = foutF$mat, FMat = foutF$FMat, LEstsM = foutF$LEstsM, ModelsM = foutF$ModelsM ) 
   }
  return(samonOut)
}
