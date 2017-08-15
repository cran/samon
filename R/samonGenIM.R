samonGenIM <- function(mat,
            Npart = 10,
            InitialSigmaH = 1.0,  
            HighSigmaH = 2.0,
            InitialSigmaF = 1.0,
            HighSigmaF = 2.0,
            inmodel = inmodel,
            NSamples = 0,
            seed0 = 1,
            seed1 = 1,           
            MaxIter = 25,
            FAconvg = 1E-7,
            FRconvg = 1E-7,
            SAconvg = 1E-7
                  )
{
  nt <- ncol(mat)
  n0 <- nrow(mat)
  
  # some sizes
  nmodel   <- 6
  nfills   <- 1

  FMat     <- matrix( 0, n0 * nfills, nt + 2,    )
  LEstsM   <- matrix( 0, nt * nfills, nmodel +  7)
  ModelsM  <- matrix( 0, nt * nfills, nmodel     )

  HM       <- matrix( 0, nrow =   nfills,  ncol =    6 )
  FM       <- matrix( 0, nrow =   nfills,  ncol =    6 )
  Sample   <- matrix( 0, nrow = NSamples,  ncol = nt+1 )

  mat      <- as.matrix(mat)
  
  fout <-  .C("samon_genIM",
  n0                = as.integer (                       n0   ),
  nt                = as.integer (                       nt   ),  
  mat               = as.double  ( as.vector(           mat ) ),
              
  nmodel            = as.integer (                   nmodel   ),
  inmodel           = as.integer ( as.vector(       inmodel ) ),

  FMat              = as.double  ( as.vector(          FMat ) ),
  LEstsM            = as.double  ( as.vector(        LEstsM ) ),
  ModelsM           = as.integer ( as.vector(       ModelsM ) ),
              
  HM                = as.double  ( as.vector(            HM ) ),
  FM                = as.double  ( as.vector(            FM ) ),
              
  seed0             = as.integer (                    seed0   ),
  seed1             = as.integer (                    seed1   ),
              
  InitialSigmaH     = as.double  (            InitialSigmaH   ),  
  HighSigmaH        = as.double  (               HighSigmaH   ),
  InitialSigmaF     = as.double  (            InitialSigmaF   ),
  HighSigmaF        = as.double  (               HighSigmaF   ),
              
  Npart             = as.integer (                    Npart   ),
  NSamples          = as.integer (                 NSamples   ),

  MaxIter           = as.integer (                  MaxIter   ),
  FAconvg           = as.double  (                  FAconvg   ),
  FRconvg           = as.double  (                  FRconvg   ),
  SAconvg           = as.double  (                  SAconvg   ),
              
  NFills            = as.integer (                   nfills   ),
              
  Sample            = as.double  ( as.vector(        Sample ) ),
  NAOK=TRUE
      )

  hnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaH", "lossH" )
  fnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaF", "lossF" )
  anames          <- c("Sample", "Type", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
  IFinamesMB      <- c("Sample", "alpha", "Obsno", "AEst", "IFEst" )

  # always done 
  if ( TRUE ) {
    fout$HM              <- matrix( fout$HM,     nrow = nfills,     ncol =        6,  byrow = TRUE  )
    fout$FM              <- matrix( fout$FM,     nrow = nfills,     ncol =        6,  byrow = TRUE  )
    colnames(fout$HM)    <- hnames
    colnames(fout$FM)    <- fnames
    
    fout$inmodel         <- matrix(fout$inmodel, nrow= nt,          ncol =   nmodel,  byrow = FALSE )
    
    fout$FMat            <- matrix(fout$FMat,    nrow= n0 * nfills, ncol =   nt + 2,  byrow = TRUE  )
    fout$LEstsM          <- matrix(fout$LEstsM,  nrow= nt * nfills, ncol = nmodel+7,  byrow = TRUE  )
    fout$ModelsM         <- matrix(fout$ModelsM, nrow= nt * nfills, ncol =   nmodel,  byrow = TRUE  )
  }
  
  if ( NSamples > 0 ) {
    fout$Sample          <- matrix(fout$Sample,  nrow=    NSamples, ncol =      nt+1, byrow = TRUE  )
    ## change NaN to NA
    fout$Sample[ is.nan(fout$Sample) ] <- NA
  }
  
  fout$mat <- matrix(fout$mat, nrow= n0,ncol = nt, byrow=FALSE  )
  return(fout)
}
