samonIM <- function(mat, Npart = 10, InitialSigmaH = 1.0, HighSigmaH = 2.0, InitialSigmaF = 1.0, HighSigmaF = 2.0,
           inmodel = inmodel, NSamples = 0, NIMimpute = 1, lb = 0, ub = 101, zeta1 = 1, zeta2 = 1, seed0 = 1,
           seed1 = 1, MaxIter = 25, FAconvg = 1E-7, FRconvg = 1E-7, SAconvg = 1E-7, alphaList = c(0), retIFiM  = FALSE,
           retIFiS  = FALSE, retSample = FALSE, retFMatM = FALSE, retFMatS = FALSE, Tfun = NULL      
           )
{
    
  nt         <- ncol(mat)
  n0         <- nrow(mat)

  if ( !is.null(Tfun) ) lenTfun <- nrow(Tfun)
  else lenTfun <- 0
  
  # convert to numeric indicators 0/1
  retIFiM    <- as.numeric(retIFiM)
  retIFiS    <- as.numeric(retIFiS)
  retFMatM   <- as.numeric(retFMatM)
  retFMatS   <- as.numeric(retFMatS)
  retSample  <- as.numeric(retSample)
  
  SampleSize <- NSamples * n0
  nfills     <- NIMimpute + 1
  nalpha     <- length( alphaList )
  
  # compute some row sizes 
  nrowsM       <- c( 1, 1, nalpha ) * nfills
  nrowsS       <- c( 1, 1, nalpha )
  nrowsIFiM    <- 1
  nrowsIFiS    <- 1
  nrowsFMatM   <- 1
  nrowsFMatS   <- 1

  if ( retFMatM == 1 ) nrowsFMatM <- n0 * nfills
  if ( retFMatS == 1 ) nrowsFMatS <- n0 * nfills * NSamples
  
  # max model size
  nmodel     <- 6
  FMatM      <- matrix(0,  nrowsFMatM,     nt + 2)
  LEstsM     <- matrix(0, nt * nfills, nmodel + 7)
  ModelsM    <- matrix(0, nt * nfills, nmodel    )

  FMatS      <- matrix(0,             nrowsFMatS,     nt + 2)
  LEstsS     <- matrix(0, nt * nfills * NSamples, nmodel + 7)
  ModelsS    <- matrix(0, nt * nfills * NSamples, nmodel    )

  if (retIFiM   == 1 ) nrowsIFiM   <- n0 * nalpha * nfills
  if (retIFiS   == 1 ) nrowsIFiS   <- n0 * nalpha * nfills * NSamples
  if (NSamples  >  0 ) nrowsS      <- nrowsS * nfills * NSamples

  HM      <- matrix(0,  nrow = nrowsM[1],  ncol =  6 )
  FM      <- matrix(0,  nrow = nrowsM[2],  ncol =  6 )
  IFM     <- matrix(0,  nrow = nrowsM[3],  ncol =  7 )

  IFiM    <- matrix(0,  nrow = nrowsIFiM,  ncol =  5 )

  HS      <- matrix(0,  nrow = nrowsS[1],  ncol =  6 )
  FS      <- matrix(0,  nrow = nrowsS[2],  ncol =  6 )
  IFS     <- matrix(0,  nrow = nrowsS[3],  ncol =  7 )

  IFiS    <- matrix(0,  nrow = nrowsIFiS,  ncol =  5 )
  
  if ( retSample == 1 ) Sample  <- matrix(0,  nrow = SampleSize, ncol = nt + 1)
  else Sample <- matrix(0, nrow = 0, ncol = nt + 1)
  
  mat     <- as.matrix(mat)
  
  fout <-  .C("samon_ngenIMIF",
  n0                = as.integer (                       n0   ),
  nt                = as.integer (                       nt   ),  
  mat               = as.double  ( as.vector(           mat ) ),
              
  nmodel            = as.integer (                   nmodel   ),
  inmodel           = as.integer ( as.vector(       inmodel ) ),
              
  FMatM             = as.double  ( as.vector(         FMatM ) ),
  LEstsM            = as.double  ( as.vector(        LEstsM ) ),
  ModelsM           = as.integer ( as.vector(       ModelsM ) ),
              
  FMatS             = as.double  ( as.vector(         FMatS ) ),
  LEstsS            = as.double  ( as.vector(        LEstsS ) ),
  ModelsS           = as.integer ( as.vector(       ModelsS ) ),
              
  HM                = as.double  ( as.vector(            HM ) ),
  FM                = as.double  ( as.vector(            FM ) ),
  IFM               = as.double  ( as.vector(           IFM ) ),
              
  HS                = as.double  ( as.vector(            HS ) ),
  FS                = as.double  ( as.vector(            FS ) ),
  IFS               = as.double  ( as.vector(           IFS ) ),

  seed0             = as.integer (                    seed0   ),
  seed1             = as.integer (                    seed1   ),
              
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

  Nfills            = as.integer (                   nfills   ),
              
  retIFiM           = as.integer (                  retIFiM   ),
  retIFiS           = as.integer (                  retIFiS   ),
  IFiM              = as.double  ( as.vector(          IFiM ) ),
  IFiS              = as.double  ( as.vector(          IFiS ) ),

  retSample         = as.integer (                retSample   ),
  Sample            = as.double  ( as.vector(        Sample ) ),

  retFMatM          = as.integer (                 retFMatM   ),
  retFMatS          = as.integer (                 retFMatS   ),

  lenTfun           = as.integer (                  lenTfun   ),
  Tfun              = as.double  ( as.vector(          Tfun ) ),            
  NAOK=TRUE
      )

  hnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaH", "lossH" )
  fnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaF", "lossF" )
  anames          <- c("Sample", "Type", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
  anamesjk        <- c("Sample", "Dropped", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
  IFinamesMB      <- c("Sample", "alpha", "Obsno", "AEst", "IFEst" )

  if ( TRUE ) {
    fout$HM              <- matrix( fout$HM,     nrow = nrowsM[1],  ncol =  6,       byrow = TRUE )
    fout$FM              <- matrix( fout$FM,     nrow = nrowsM[2],  ncol =  6,       byrow = TRUE )
    fout$IFM             <- matrix( fout$IFM,    nrow = nrowsM[3],  ncol =  7,       byrow = TRUE )
    colnames(fout$HM)    <- hnames
    colnames(fout$FM)    <- fnames
    colnames(fout$IFM)   <- anames

    if ( retFMatM == 1 ) fout$FMatM <- matrix(fout$FMatM, nrow = nrowsFMatM, ncol = nt + 2, byrow=TRUE   )
    else fout$FMatM <- NULL
    
    fout$LEstsM          <- matrix(fout$LEstsM,  nrow= nt * nfills, ncol = nmodel+7,  byrow=TRUE   )
    fout$ModelsM         <- matrix(fout$ModelsM, nrow= nt * nfills, ncol = nmodel,    byrow=TRUE   )
  }
  
  if ( NSamples > 0 ) {
    if ( retSample == 1 ) fout$Sample <- matrix( fout$Sample, nrow = SampleSize, ncol = nt + 1,   byrow =  TRUE  )
    else fout$Sample <- NULL
    
    fout$HS              <- matrix( fout$HS,     nrow = nrowsS[1],  ncol =  6,       byrow = TRUE )
    fout$FS              <- matrix( fout$FS,     nrow = nrowsS[2],  ncol =  6,       byrow = TRUE )
    fout$IFS             <- matrix( fout$IFS,    nrow = nrowsS[3],  ncol =  7,       byrow = TRUE )
    colnames(fout$HS)    <- hnames
    colnames(fout$FS)    <- fnames
    colnames(fout$IFS)   <- anames

    if ( retFMatS == 1 ) fout$FMatS <- matrix(fout$FMatS, nrow = nrowsFMatS, ncol = nt + 2, byrow=TRUE   )
    else fout$FMatS <- NULL
    fout$LEstsS          <- matrix(fout$LEstsS,  nrow= nt * NSamples * nfills,  ncol = nmodel+7, byrow=TRUE   )
    fout$ModelsS         <- matrix(fout$ModelsS, nrow= nt * NSamples * nfills,  ncol = nmodel,   byrow=TRUE   )
    
  } else {
    fout$Sample <- NULL
    fout$HS     <- fout$FS <- fout$IFS <- fout$IFMatS <- fout$LEstsS <- fout$IFS <- NULL
  }
  fout$retFMatM   <- fout$retFMatS <- fout$retSample <- fout$retIFiM <- fout$retIFiS <- NULL 
  
  if ( lenTfun > 0 ) {
     fout$Tfun           <- matrix( fout$Tfun, nrow = lenTfun, ncol = 2, byrow = FALSE )
     colnames(fout$Tfun) <- c("uniqueValue", "TiltValue")
  } else {
     fout$Tfun <- NULL
  }
  fout$NIMimpute <- NIMimpute
  fout$mat       <- matrix(fout$mat, nrow= n0, ncol = nt, byrow=FALSE  )
  fout$inmodel   <- inmodel
  return(fout)
}
