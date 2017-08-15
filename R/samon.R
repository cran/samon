samon <- function(mat,
            Npart = 10,
            InitialSigmaH = 1.0,  
            HighSigmaH = 2.0,
            InitialSigmaF = 1.0,
            HighSigmaF = 2.0,
            lb = 0,
            ub = 101,
            zeta1 = 1,
            zeta2 = 1,
            NSamples = 0,
            seed0 = 1,
            MaxIter = 25,
            FAconvg = 1E-7,
            FRconvg = 1E-7,
            SAconvg = 1E-7,
            alphaList = c(0),
            MJackknife = FALSE,
            SJackknife = FALSE,      
            retIFiM  = FALSE,
            retIFiS  = FALSE,
            Tfun     = NULL      
                  )
{
    
  nt <- ncol(mat)
  n0 <- nrow(mat)
  nv <- nt*n0


  if ( !is.null(Tfun) ) {
      lenTfun <- nrow(Tfun)
  } else {
      lenTfun <- 0
  }
  nalpha <- length( alphaList )
  
  # convert to numeric indicators 0/1
  MJackknife <- as.numeric(MJackknife)
  SJackknife <- as.numeric(SJackknife)
  retIFiM    <- as.numeric(retIFiM)
  retIFiS    <- as.numeric(retIFiS)

  # some imply the others
  if ( NSamples   <= 0 ) SJackknife <- 0

  # compute some row sizes 
  nrowsM       <- c( 1, 1, nalpha )
  nrowsDefault <- c( 1, 1,      1 )
  nrowsMjk     <- nrowsDefault
  nrowsS       <- nrowsDefault
  nrowsSjk     <- nrowsDefault
  nrowsIFiM    <- 1
  nrowsIFiS    <- 1
  
  if (MJackknife == 1 ) nrowsMjk    <- nrowsM * n0
  if (NSamples    > 0 ) nrowsS      <- nrowsM * NSamples
  if (SJackknife == 1 ) nrowsSjk    <- nrowsM * n0 * NSamples
  if (retIFiM    == 1 ) nrowsIFiM   <- n0 * nalpha
  if (retIFiS    == 1 ) nrowsIFiS   <- n0 * nalpha * NSamples

  HM      <- matrix(0,  nrow = nrowsM[1],      ncol =  6 )
  FM      <- matrix(0,  nrow = nrowsM[2],      ncol =  6 )
  IFM     <- matrix(0,  nrow = nrowsM[3],      ncol =  7 )

  HMjk    <- matrix(0,  nrow = nrowsMjk[1],    ncol =  6 )
  FMjk    <- matrix(0,  nrow = nrowsMjk[2],    ncol =  6 )
  IFMjk   <- matrix(0,  nrow = nrowsMjk[3],    ncol =  7 )

  HS      <- matrix(0,  nrow = nrowsS[1],      ncol =  6 )
  FS      <- matrix(0,  nrow = nrowsS[2],      ncol =  6 )
  IFS     <- matrix(0,  nrow = nrowsS[3],      ncol =  7 )
  
  HSjk    <- matrix(0,  nrow = nrowsSjk[1],    ncol =  6 )
  FSjk    <- matrix(0,  nrow = nrowsSjk[2],    ncol =  6 )
  IFSjk   <- matrix(0,  nrow = nrowsSjk[3],    ncol =  7 )

  IFiM    <- matrix(0,  nrow = nrowsIFiM,      ncol =  5 )
  IFiS    <- matrix(0,  nrow = nrowsIFiS,      ncol =  5 )

  mat     <- as.matrix(mat)
  
  fout <-  .C("samon_boot_jk2",
  n0                = as.integer (                       n0   ),
  nt                = as.integer (                       nt   ),  
  mat               = as.double  ( as.vector(           mat ) ),
              
  HM                = as.double  ( as.vector(            HM ) ),
  FM                = as.double  ( as.vector(            FM ) ),
  IFM               = as.double  ( as.vector(           IFM ) ),

  HMjk              = as.double  ( as.vector(          HMjk ) ),
  FMjk              = as.double  ( as.vector(          FMjk ) ),
  IFMjk             = as.double  ( as.vector(         IFMjk ) ),

  HS                = as.double  ( as.vector(            HS ) ),
  FS                = as.double  ( as.vector(            FS ) ),
  IFS               = as.double  ( as.vector(           IFS ) ),

  HSjk              = as.double  ( as.vector(          HSjk ) ),
  FSjk              = as.double  ( as.vector(          FSjk ) ),
  IFSjk             = as.double  ( as.vector(         IFSjk ) ),
              
  seed0             = as.integer (                    seed0   ),
              
  InitialSigmaH     = as.double  (            InitialSigmaH   ),  
  HighSigmaH        = as.double  (               HighSigmaH   ),
  InitialSigmaF     = as.double  (            InitialSigmaF   ),
  HighSigmaF        = as.double  (               HighSigmaF   ),
              
  lb                = as.double  (                       lb   ),
  ub                = as.double  (                       ub   ),
  zeta1             = as.double  (                    zeta1   ),
  zeta2             = as.double  (                    zeta2   ),
              
  Npart             = as.integer (                    Npart   ),
  NSamples          = as.integer (                 NSamples   ),

  MaxIter           = as.integer (                  MaxIter   ),
  FAconvg           = as.double  (                  FAconvg   ),
  FRconvg           = as.double  (                  FRconvg   ),
  SAconvg           = as.double  (                  SAconvg   ),
              
  Nalpha            = as.integer (                   nalpha   ),
  alphaList         = as.double  (                alphaList   ),
              
  MJackknife        = as.integer (               MJackknife   ),
  SJackknife        = as.integer (               SJackknife   ),
              
  retIFiM           = as.integer (                  retIFiM   ),
  retIFiS           = as.integer (                  retIFiS   ),

  IFiM              = as.double  ( as.vector(          IFiM ) ),
  IFiS              = as.double  ( as.vector(          IFiS ) ),

  lenTfun           = as.integer (                  lenTfun   ),
  Tfun              = as.double  ( as.vector(          Tfun ) ),            
  NAOK=TRUE
      )

  hnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaH", "lossH" )
  fnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaF", "lossF" )
  anames          <- c("Sample", "Type", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
  anamesjk        <- c("Sample", "Dropped", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
  IFinamesMB      <- c("Sample", "alpha", "Obsno", "AEst", "IFEst" )

  # always done 
  if ( TRUE ) {
    fout$HM              <- matrix( fout$HM,     nrow = nrowsM[1],     ncol =  6, byrow = TRUE )
    fout$FM              <- matrix( fout$FM,     nrow = nrowsM[2],     ncol =  6, byrow = TRUE )
    fout$IFM             <- matrix( fout$IFM,    nrow = nrowsM[3],     ncol =  7, byrow = TRUE )
    colnames(fout$HM)    <- hnames
    colnames(fout$FM)    <- fnames
    colnames(fout$IFM)   <- anames
  }
  
  if ( NSamples > 0 ) {
    fout$HS              <- matrix( fout$HS,     nrow = nrowsS[1],     ncol =  6, byrow = TRUE )
    fout$FS              <- matrix( fout$FS,     nrow = nrowsS[2],     ncol =  6, byrow = TRUE )
    fout$IFS             <- matrix( fout$IFS,    nrow = nrowsS[3],     ncol =  7, byrow = TRUE )
    colnames(fout$HS)    <- hnames
    colnames(fout$FS)    <- fnames
    colnames(fout$IFS)   <- anames
  } else {
    fout$HS              <- NULL
    fout$FS              <- NULL
    fout$IFS             <- NULL
  }

  if ( MJackknife == 1 ) {
    fout$HMjk            <- matrix( fout$HMjk,   nrow = nrowsMjk[1],   ncol =  6, byrow = TRUE )
    fout$FMjk            <- matrix( fout$FMjk,   nrow = nrowsMjk[2],   ncol =  6, byrow = TRUE )
    fout$IFMjk           <- matrix( fout$IFMjk,  nrow = nrowsMjk[3],   ncol =  7, byrow = TRUE )
    colnames(fout$HMjk)  <- hnames
    colnames(fout$FMjk)  <- fnames
    colnames(fout$IFMjk) <- anamesjk
  } else {
    fout$HMjk            <- NULL
    fout$FMjk            <- NULL
    fout$IFMjk           <- NULL
  }

  if ( SJackknife == 1 ) {
    fout$HSjk            <- matrix( fout$HSjk,   nrow = nrowsSjk[1],   ncol =  6, byrow = TRUE )
    fout$FSjk            <- matrix( fout$FSjk,   nrow = nrowsSjk[2],   ncol =  6, byrow = TRUE )
    fout$IFSjk           <- matrix( fout$IFSjk,  nrow = nrowsSjk[3],   ncol =  7, byrow = TRUE )
    colnames(fout$HSjk)  <- hnames
    colnames(fout$FSjk)  <- fnames
    colnames(fout$IFSjk) <- anamesjk
  } else {
    fout$HSjk            <- NULL
    fout$FSjk            <- NULL
    fout$IFSjk           <- NULL
  }
  
  if ( retIFiM == 1 ) {
    fout$IFiM <- matrix( fout$IFiM,   nrow = nrowsIFiM,    ncol =  5, byrow = TRUE )
    colnames(fout$IFiM)  <- IFinamesMB
  } else {
    fout$IFiM <- NULL
  }
  if ( retIFiS == 1 ) {
    fout$IFiS <- matrix( fout$IFiS,   nrow = nrowsIFiS,    ncol =  5, byrow = TRUE )
    colnames(fout$IFiS)  <- IFinamesMB
  } else {
    fout$IFiS <- NULL
  }

  if ( lenTfun > 0 ) {
     fout$Tfun           <- matrix( fout$Tfun, nrow = lenTfun, ncol = 2, byrow = FALSE )
     colnames(fout$Tfun) <- c("uniqueValue", "TiltValue")
  } else {
     fout$Tfun <- NULL
  }
  
  fout$mat        <- matrix(fout$mat,     nrow= n0,          ncol = nt, byrow=FALSE  )
    
  return(fout)
}
