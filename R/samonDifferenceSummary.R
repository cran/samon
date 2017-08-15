# --------------------------------------------------------------------------------
# Summary function for samon difference in treatment
# --------------------------------------------------------------------------------
samonDifferenceSummary <- function(trt1,trt2,CIlevel=0.95) {

  TM1 <- trt1$TM[, c("alpha","IFEst","IFVar","meanjk","sejk","meanBoot","seBoot")]
  TM2 <- trt2$TM[, c("alpha","IFEst","IFVar","meanjk","sejk","meanBoot","seBoot")]
  TS1 <- trt1$TS[, c("alpha","Sample","IFEst","IFVar","meanjk","sejk","MIFEst")]
  TS2 <- trt2$TS[, c("alpha","Sample","IFEst","IFVar","meanjk","sejk","MIFEst")]

  Nalpha    <- trt1$Nalpha
  alphaList <- trt1$alphaList
  NSamples  <- trt1$NSamples
  n10       <- trt1$n0
  n20       <- trt2$n0

  colnames(TM1) <- c("alpha","IFEst1","IFVar1","meanjk1","sejk1","meanBoot1","seBoot1")
  colnames(TM2) <- c("alpha","IFEst2","IFVar2","meanjk2","sejk2","meanBoot2","seBoot2")
  colnames(TS1) <- c("alpha","Sample","IFEst1","IFVar1","meanjk1","sejk1","MIFEst1")
  colnames(TS2) <- c("alpha","Sample","IFEst2","IFVar2","meanjk2","sejk2","MIFEst2")

  TM                <- merge(TM1,TM2,by=c("alpha"))
  TM[,"Difference"] <- TM[,"IFEst2" ] - TM[,"IFEst1" ]
  TM[,"DIFSE"]      <- sqrt(TM[,"IFVar2"] + TM[,"IFVar1"])
  TM[,"DjkSE"]      <- sqrt(TM[,"sejk2"]^2 + TM[,"sejk1"]^2)
  
  TS <- merge(TS1,TS2,by=c("Sample","alpha"))
  TS[,"Difference"] <- TS[,"IFEst2" ] - TS[,"IFEst1" ]
  TS[,"DIFSE"]      <- sqrt(TS[,"IFVar2"] + TS[,"IFVar1"])
  TS[,"DjkSE"]      <- sqrt(TS[,"sejk2"]^2 + TS[,"sejk1"]^2)
  TS[,"TIF"]        <- (TS[,"Difference"] - TS[,"MIFEst2"] + TS[,"MIFEst1"]) / TS[,"DIFSE"]
  TS[,"TVals"]      <- (TS[,"Difference"] - TS[,"MIFEst2"] + TS[,"MIFEst1"]) / TS[,"DjkSE" ]
  
  TS    <- data.frame(TS)
  TS[,"alpha"] <- factor(TS[,"alpha"])

  ## quantiles
  cicut  <- CIlevel
  cicut1 <- (1 - CIlevel)/2
  cicut2 <- 1 - cicut1
  myq <- function(x) {
    x1 <- as.vector(quantile(x, c(cicut1,cicut2)))
    x2 <- as.vector(quantile(abs(x), c(cicut)))
    return(c(x1,x2))                
  }
  ag1 <- aggregate( . ~ alpha, TS[,c("alpha","TVals")],      myq )
  ag2 <- aggregate( . ~ alpha, TS[,c("alpha","Difference")], myq )
  ag3 <- aggregate( . ~ alpha, TS[,c("alpha","TIF")  ],      myq )

  ag1 <- cbind(ag1[1],ag1[[2]])  
  ag2 <- cbind(ag2[1],ag2[[2]])  
  ag3 <- cbind(ag3[1],ag3[[2]])
  
  colnames(ag1) <- c("alpha", "tLow",   "tHigh",     "tSym")
  colnames(ag2) <- c("alpha", "DIFLow", "DIFHigh",   "DIFSym")
  colnames(ag3) <- c("alpha", "tIFLow", "tIFHigh",   "tIFSym")

  TM <- merge(TM,ag1,by=c("alpha"))
  TM <- merge(TM,ag2,by=c("alpha"))
  TM <- merge(TM,ag3,by=c("alpha"))

  zcut <- qnorm( cicut2 )
  ## IF CI
  lb1 <- TM[,"Difference"] - zcut * TM[,"DIFSE"]
  ub1 <- TM[,"Difference"] + zcut * TM[,"DIFSE"]

  ## Bootstrap CI
  lb2 <- TM[,"Difference"] - zcut * TM[,"DIFSE" ]
  ub2 <- TM[,"Difference"] + zcut * TM[,"DIFSE" ]

  ## jk CI
  lb3 <- TM[,"Difference"] - zcut * TM[,"DjkSE" ]
  ub3 <- TM[,"Difference"] + zcut * TM[,"DjkSE" ]

  ## quantiles
  lb4 <- TM[,"DIFLow" ]
  ub4 <- TM[,"DIFHigh"]

  ## t using IF se
  lb5 <- TM[,"Difference"] - TM[,"tIFHigh"] * TM[,"DIFSE"]
  ub5 <- TM[,"Difference"] - TM[,"tIFLow" ] * TM[,"DIFSE"]

  ## t using jk se
  lb6 <- TM[,"Difference"] - TM[,"tHigh"] * TM[,"DjkSE"]
  ub6 <- TM[,"Difference"] - TM[,"tLow" ] * TM[,"DjkSE"]

  ## t using IF se symmetric
  lb7 <- TM[,"Difference"] - TM[,"tIFSym"] * TM[,"DIFSE"]
  ub7 <- TM[,"Difference"] + TM[,"tIFSym"] * TM[,"DIFSE"]
  
  ## t using jk se symmetric
  lb8 <- TM[,"Difference"] - TM[,"tIFSym"] * TM[,"DjkSE"]
  ub8 <- TM[,"Difference"] + TM[,"tIFSym"] * TM[,"DjkSE"]

  CI <- cbind( alphaList, lb1, ub1, lb2, ub2, lb3, ub3, lb4, ub4, lb5, ub5, lb6, ub6, lb7, ub7, lb8, ub8)
  
  Ret <- list( TM = TM, CI = CI, n10 = n10, n20 = n20, NSamples = NSamples, Nalpha = Nalpha, alphaList = alphaList, CIlevel=CIlevel )
  return(Ret)
  
}
