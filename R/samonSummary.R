# --------------------------------------------------------------------------------
# Summary function for samon objects
# --------------------------------------------------------------------------------
samonSummary <- function(trt, CIlevel=0.95) {

  trtnames <- names(trt)
  
  HM       <- trt$HM   [,c("Convergence","Iterations","SigmaH","lossH")]
  FM       <- trt$FM   [,c("Convergence","Iterations","SigmaF","lossF")]
  TM       <- trt$IFM  [,c("alpha","AEst","AVar","IFEst","IFVar")]  
  TMjk     <- trt$IFMjk[,c("alpha","IFEst","IFVar")]
  
  TS       <- trt$IFS  [,c("Sample","alpha","IFEst","IFVar")]
  if ( match("RegSjk", trtnames, -1 ) > 0 ) {
      TSjk <- trt$RegSjk
      colnames(TSjk) <- c("Sample","alpha","IFEst","IFVar")
      RegSjkInd      <- 1
  } else {
      TSjk  <- trt$IFSjk[,c("Sample","alpha","IFEst","IFVar")]
      RegSjkInd      <- 0
  }
  
  Nalpha    <- trt$Nalpha
  alphaList <- trt$alphaList
  NSamples  <- trt$NSamples
  n0        <- trt$n0

  TS    <- data.frame(TS)
  TS[,"alpha"] <- factor(TS[,"alpha"])

  TMjk  <- data.frame(TMjk)
  TMjk[, "alpha"] <- factor(TMjk[,"alpha"])

  ## add jk mean and se to main
  xy <- aggregate( . ~ alpha, TMjk[,c("alpha","IFEst")], function(x) { c(mean(x),var(x)) } )
  xy <- cbind(xy[1],xy[[2]])
  colnames(xy) <- c("alpha","meanjk","sejk")
  xy[,"sejk"] <- sqrt((n0-1)*(n0-1)*xy[,"sejk"]/n0)

  TM <- merge(TM,xy,by="alpha")

  ## add bootstrap mean and se to main
  xy <- aggregate( . ~ alpha, TS[,c("alpha","IFEst")], function(x) { c(mean(x),var(x)) } )
  xy <- cbind(xy[1],xy[[2]])
  colnames(xy) <- c("alpha","meanBoot","seBoot")
  xy[,"seBoot"] <- sqrt((n0-1)*xy[,"seBoot"]/n0)

  TM <- merge(TM,xy,by="alpha")

  ## mean and se from jackknife for bootstrap results and add them to bootstraps
  if ( RegSjkInd == 0 ) {
      TSjk  <- data.frame(TSjk)
      TSjk[, "Sample" ] <- factor(TSjk[,"Sample" ])
      TSjk[, "alpha"  ] <- factor(TSjk[,"alpha"  ])

      xz <- aggregate( . ~ Sample + alpha, TSjk[,c("Sample","alpha","IFEst")], function(x) { c(mean(x),var(x)) } )
      xz <- cbind(xz[1],xz[2],xz[[3]])
      colnames(xz) <- c("Sample","alpha","meanjk","sejk" )
      xz[,"sejk"] <- sqrt((n0-1)*(n0-1)*xz[,"sejk"]/n0)
  } else {
      xz <- TSjk
      colnames(xz) <- c("Sample","alpha","meanjk","sejk" )
      xz[,"sejk"] <- sqrt((n0-1)*xz[,"sejk"])
  }

  TS <- merge(TS,xz,by=c("Sample","alpha"))

  ## add main means to bootstraps
  TMSelect <- TM[,c("alpha","IFEst")]
  colnames(TMSelect) <- c("alpha","MIFEst")

  TS <- merge(TS,TMSelect,by=c("alpha"))

  ## calculate a couple of t values
  TS[,"Tvals"] <- ( TS[,"IFEst"] - TS[,"MIFEst"] ) / TS[,"sejk"  ]
  TS[,"TIF"  ] <- ( TS[,"IFEst"] - TS[,"MIFEst"] ) / sqrt(TS[,"IFVar"])

  ## quantiles
  cicut  <- CIlevel
  cicut1 <- (1 - CIlevel)/2
  cicut2 <- 1 - cicut1
  myq <- function(x) {
    x1 <- as.vector(quantile(x, c(cicut1, cicut2)))
    x2 <- as.vector(quantile(abs(x), c(cicut)))
    return(c(x1,x2))                
  }
  ag1 <- aggregate( . ~ alpha, TS[,c("alpha","Tvals")], myq )
  ag2 <- aggregate( . ~ alpha, TS[,c("alpha","IFEst")], myq )
  ag3 <- aggregate( . ~ alpha, TS[,c("alpha","TIF")  ], myq )

  ag1 <- cbind(ag1[1],ag1[[2]])  
  ag2 <- cbind(ag2[1],ag2[[2]])  
  ag3 <- cbind(ag3[1],ag3[[2]])
  
  colnames(ag1) <- c("alpha", "tLow",   "tHigh",   "tSym")
  colnames(ag2) <- c("alpha", "IFLow",  "IFHigh",  "IFSym")
  colnames(ag3) <- c("alpha", "tIFLow", "tIFHigh", "tIFSym")

  TM <- merge(TM,ag1,by=c("alpha"))
  TM <- merge(TM,ag2,by=c("alpha"))
  TM <- merge(TM,ag3,by=c("alpha"))

  zcut <- qnorm( cicut2 )
  ## IF CI
  lb1 <- TM[,"IFEst"] - zcut * sqrt(TM[,"IFVar"])
  ub1 <- TM[,"IFEst"] + zcut * sqrt(TM[,"IFVar"])

  ## using bootstrap se
  lb2 <- TM[,"IFEst"] - zcut * TM[,"seBoot"]
  ub2 <- TM[,"IFEst"] + zcut * TM[,"seBoot"]

  ## using jk se
  lb3 <- TM[,"IFEst"] - zcut * TM[,"sejk"]
  ub3 <- TM[,"IFEst"] + zcut * TM[,"sejk"]

  ## Bootstrap CI
  lb4 <- TM[,"IFLow"]
  ub4 <- TM[,"IFHigh"]

  ## t using IF SE
  lb5 <- TM[,"IFEst"] - TM[,"tIFHigh"] * sqrt(TM[,"IFVar"])
  ub5 <- TM[,"IFEst"] - TM[,"tIFLow" ] * sqrt(TM[,"IFVar"])

  ## t using jk se
  lb6 <- TM[,"IFEst"] - TM[,"tHigh"] * TM[,"sejk"]
  ub6 <- TM[,"IFEst"] - TM[,"tLow" ] * TM[,"sejk"]

  ## t using IF se symmetric
  lb7 <- TM[,"IFEst"] - TM[,"tIFSym"] * sqrt(TM[,"IFVar"])
  ub7 <- TM[,"IFEst"] + TM[,"tIFSym"] * sqrt(TM[,"IFVar"])

  ## t using jk se Symmetric
  lb8 <- TM[,"IFEst"] - TM[,"tSym"] * TM[,"sejk"]
  ub8 <- TM[,"IFEst"] + TM[,"tSym"] * TM[,"sejk"]

  CI <- cbind(alphaList, lb1, ub1, lb2, ub2, lb3, ub3, lb4, ub4, lb5, ub5, lb6, ub6, lb7, ub7, lb8, ub8)
  
  Ret <- list( TM = TM, TS = TS, CI = CI, n0 = n0, NSamples = NSamples, Nalpha = Nalpha, alphaList = alphaList, HM = HM, FM = FM, CIlevel=CIlevel )
  return(Ret)
}
