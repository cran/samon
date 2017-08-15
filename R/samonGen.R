samonGen <- function(mat,
            Npart = 10,
            InitialSigmaH = 1.0,  
            HighSigmaH = 2.0,
            InitialSigmaF = 1.0,
            HighSigmaF = 2.0,
            NSamples = 0,
            seed0 = 1,
            MaxIter = 25,
            FAconvg = 1E-7,
            FRconvg = 1E-7,
            SAconvg = 1E-7
                  )
{
  nt <- ncol(mat)
  n0 <- nrow(mat)
  
  HM      <- matrix(0,  nrow =        1,   ncol =  6 )
  FM      <- matrix(0,  nrow =        1,   ncol =  6 )
  Sample  =  matrix(0,  nrow = NSamples,   ncol = nt )

  mat     <- as.matrix(mat)
  
  fout <-  .C("samon_gen",
  n0                = as.integer (                       n0   ),
  nt                = as.integer (                       nt   ),  
  mat               = as.double  ( as.vector(           mat ) ),
              
  HM                = as.double  ( as.vector(            HM ) ),
  FM                = as.double  ( as.vector(            FM ) ),

  seed0             = as.integer (                    seed0   ),
              
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
              
  Sample            = as.double  ( as.vector(        Sample ) ),
  NAOK=TRUE
      )

    hnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaH", "lossH" )
    fnames          <- c("Sample", "Type", "Convergence", "Iterations", "SigmaF", "lossF" )
##  anames          <- c("Sample", "Type", "alpha", "AEst", "AVar", "IFEst", "IFVar" )
##  IFinamesMB      <- c("Sample", "alpha", "Obsno", "AEst", "IFEst" )

  # always done 
  if ( TRUE ) {
    fout$HM              <- matrix( fout$HM,  nrow = 1,   ncol =  6, byrow = TRUE )
    fout$FM              <- matrix( fout$FM,  nrow = 1,   ncol =  6, byrow = TRUE )
    colnames(fout$HM)    <- hnames
    colnames(fout$FM)    <- fnames
  }
  
  if ( NSamples > 0 ) {
    fout$Sample          <- matrix( fout$Sample,  nrow = NSamples,     ncol = nt, byrow = TRUE )
    ## change NaN to NA
    fout$Sample[ is.nan(fout$Sample) ] <- NA
  }
  
  fout$mat <- matrix(fout$mat, nrow= n0, ncol = nt, byrow=FALSE  )
    
  return(fout)
}
