alpha0TableResultsIM <- function( Y1, samonR1, Y2, samonR2, samonRD ) {

  # original data at last time-point and completers
  Y1L <- Y1[, ncol(Y1)]
  Y2L <- Y2[, ncol(Y2)]
  Y1C <- Y1L[ !is.na(Y1L) ]
  Y2C <- Y2L[ !is.na(Y2L) ]

  # IF results
  TM1    <- as.matrix(samonR1$TM)
  TM2    <- as.matrix(samonR2$TM)
  TMD    <- as.matrix(samonRD$TM)

  # alpha 0 TM
  TM10   <- TM1[ TM1[,1] == 0, ]
  TM20   <- TM2[ TM2[,1] == 0, ]
  TMD0   <- TMD[ TMD[,1] == 0, ]
  
  # Confidence interval results
  CI1    <- as.matrix(samonR1$CI)
  CI2    <- as.matrix(samonR2$CI)
  CID    <- as.matrix(samonRD$CI)

  # alpha 0 CI
  CI10   <- CI1[ CI1[,1] == 0, ]
  CI20   <- CI2[ CI2[,1] == 0, ]
  CID0   <- CID[ CID[,1] == 0, ]
  
  # H,F optimization results
  optP1  <- samonR1$HM[,c("Convergence", "Iterations", "SigmaH", "lossH")]
  optQ1  <- samonR1$FM[,c("Convergence", "Iterations", "SigmaF", "lossF")]
  optP2  <- samonR2$HM[,c("Convergence", "Iterations", "SigmaH", "lossH")]
  optQ2  <- samonR2$FM[,c("Convergence", "Iterations", "SigmaF", "lossF")]
  
##  optP1  <- matrix(samonR1$PM,nrow=1)
##  optQ1  <- matrix(samonR1$QM,nrow=1)
##  optP2  <- matrix(samonR2$PM,nrow=1)
##  optQ2  <- matrix(samonR2$QM,nrow=1)
  
  # Table 1, n, mean and sd at each time-point 
  Miss <- matrix(NA, 10, ncol(Y1))

  Miss[  1,]  <- apply(Y1,2, function(x) sum (  is.na(x) ))
  Miss[  2,]  <- apply(Y1,2, function(x) sum ( !is.na(x) ))
  Miss[  3,]  <- Miss[1,] + Miss[2,]
  Miss[  4,]  <- apply(Y1,2, function(x) mean( x[ !is.na(x) ]))
  Miss[  5,]  <- apply(Y1,2, function(x) sd  ( x[ !is.na(x) ]))

  Miss[  6,]  <- apply(Y2,2, function(x) sum (  is.na(x) ))
  Miss[  7,]  <- apply(Y2,2, function(x) sum ( !is.na(x) ))
  Miss[  8,]  <- Miss[6,] + Miss[7,]
  Miss[  9,]  <- apply(Y2,2, function(x) mean( x[ !is.na(x) ]))
  Miss[ 10,]  <- apply(Y2,2, function(x) sd  ( x[ !is.na(x) ]))

  Miss2 <- t(Miss)
  Miss2 <- Miss2[, c(2,4,5,7,9,10)]
  Miss2 <- cbind( 1:nrow(Miss2), Miss2 )
  colnames(Miss2) <- c("t", "Trt1N", "Trt1mean", "Trt1SD",
                            "Trt2N", "Trt2mean", "Trt2SD" )

  # tables 2 and 3, completers
  ObsRes       <- matrix(NA, 2, 6)
  ObsRes[1,]   <- c( 1, ncol(Y1), length(Y1L), length(Y1C), mean(Y1C), sd(Y1C) )
  ObsRes[2,]   <- c( 2, ncol(Y2), length(Y2L), length(Y2C), mean(Y2C), sd(Y2C) )
  colnames(ObsRes) <- c("Treatment", "NT","N0","NComplete","meanComplete","sdComplete")

  # ttest between completers
  ttest021 <- t.test(Y2C,Y1C)
  
  TtestRes21     <- matrix(NA, 1, 5)
  TtestRes21[1,] <- c( mean(Y2C) - mean(Y1C), ttest021$statistic, ttest021$conf.int[[1]], ttest021$conf.int[[2]], ttest021$p.value)
  colnames(TtestRes21) <- c("Difference", "ttestStatistic","LCL","UCL","Pvalue") 


  # tables 4 and 5 LOCF
  # last available observation carried forward
  Y1Locf <- Y1[, 1]
  Lt1    <- rep(1,nrow(Y1))
  for ( i in 1:nrow(Y1) ) {
    for ( t in 2:ncol(Y1)) {
     if ( !is.na(Y1[i,t]) ) {
       Y1Locf[i] = Y1[i,t]
       Lt1[i]    = t
     }
    }
  }

  Y2Locf <- Y2[, 1]
  Lt2    <- rep(1,nrow(Y2))
  for ( i in 1:nrow(Y2) ) {
    for ( t in 2:ncol(Y2)) {
     if ( !is.na(Y2[i,t]) ) {
       Y2Locf[i] = Y2[i,t]
       Lt2[i]    = t
     }
    }
  }

  Locf      <- matrix(NA,2,4)
  Locf[1,]  <- c( 1, sum( !is.na(Y1Locf) ), mean(Y1Locf), sd(Y1Locf))
  Locf[2,]  <- c( 2, sum( !is.na(Y2Locf) ), mean(Y2Locf), sd(Y2Locf))
  colnames(Locf) <- c("Treatment", "N","mean","sd")

  ttestLocf21 <- t.test(Y2Locf,Y1Locf)

  TtestLRes21     <- matrix(NA, 1, 5)
  TtestLRes21[1,] <- c( mean(Y2Locf) - mean(Y1Locf), ttestLocf21$statistic, ttestLocf21$conf.int[[1]], ttestLocf21$conf.int[[2]], ttestLocf21$p.value)
  colnames(TtestLRes21) <- c("Difference", "ttestStatistic","LCL","UCL","Pvalue") 

  # table 6 IF results
  IFRes        <- matrix(NA, 3, 3)
##  IFRes[1,]    <- c(TM10[  4 ], CI10[ "lb7" ], CI10[ "ub7" ], CI10[ "lb8" ], CI10[ "ub8" ])
##  IFRes[2,]    <- c(TM20[  4 ], CI20[ "lb7" ], CI20[ "ub7" ], CI20[ "lb8" ], CI20[ "ub8" ])
##  IFRes[3,]    <- c(TMD0[ 18 ], CID0[ "lb7" ], CID0[ "ub7" ], CID0[ "lb8" ], CID0[ "ub8" ])

  IFRes[1,]    <- c(TM10[      "IFEst" ], CI10[ "lb7" ], CI10[ "ub7" ] )
  IFRes[2,]    <- c(TM20[      "IFEst" ], CI20[ "lb7" ], CI20[ "ub7" ] )
  IFRes[3,]    <- c(TMD0[ "Difference" ], CID0[ "lb7" ], CID0[ "ub7" ] )
  colnames(IFRes) <- c("Estimate", "LCLBootstrap", "UCLBootstrap")

  CIlevel <- samonRD$CIlevel
  if ( is.null(CIlevel) ) CIlevel <- .95

  # table 7 Optimal sigmaP and sigmaq
  optRes <- rbind( as.matrix(cbind( 1, optP1, optQ1 )), as.matrix(cbind( 2, optP2, optQ2 )) )
  colnames(optRes) <- c("Treatment",
                        "ConvergenceH", "IterH", "SigmaH", "lossH",
                        "ConvergenceF", "IterF", "SigmaF", "lossF" )

  ReturnList <- list( Obsmeans = Miss2, ObsResults = ObsRes, ObsTtest = TtestRes21, LOCF = Locf, LOCFTtest = TtestLRes21, IFResults = IFRes, optRes = optRes, CIlevel=CIlevel )
  return(ReturnList)
}
