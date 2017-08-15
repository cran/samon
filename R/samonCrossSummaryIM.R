# --------------------------------------------------------------------------------
# Summary function for samon differences in treatment for pairs of alpha
# --------------------------------------------------------------------------------
samonCrossSummaryIM <- function(trt1,trt2, CIlevel=0.95) {

  Nalpha    <- trt1$Nalpha
  alphaList <- trt1$alphaList
  NSamples  <- trt1$NSamples
  n10       <- trt1$n0
  n20       <- trt2$n0
    
  TM1 <- trt1$TM[, c("alpha","IFEst","MIIFVar")]
  TM2 <- trt2$TM[, c("alpha","IFEst","MIIFVar")]
  TS1 <- trt1$TS[, c("alpha","Sample","IFEst","MIIFVar","mIFEst")]
  TS2 <- trt2$TS[, c("alpha","Sample","IFEst","MIIFVar","mIFEst")]

  colnames(TM1) <- c("alpha1","IFEst1","IFVar1")
  colnames(TM2) <- c("alpha2","IFEst2","IFVar2")
  colnames(TS1) <- c("alpha1","Sample","IFEst1","IFVar1","MIFEst1")
  colnames(TS2) <- c("alpha2","Sample","IFEst2","IFVar2","MIFEst2")

  for ( i1 in 1:Nalpha ) {
  for ( i2 in 1:Nalpha ) {
    alpha1 <- alphaList[i1]
    alpha2 <- alphaList[i2]

    TMA1 <- TM1[TM1[,"alpha1"] == alpha1,, drop=FALSE ]
    TMA2 <- TM2[TM2[,"alpha2"] == alpha2,, drop=FALSE ]
    TSA1 <- TS1[TS1[,"alpha1"] == alpha1,, drop=FALSE ]
    TSA2 <- TS2[TS2[,"alpha2"] == alpha2,, drop=FALSE ]
    
    TM                <- cbind(TMA1,TMA2)
    TM[,"Difference"] <- TM[,"IFEst2" ] - TM[,"IFEst1" ]
    TM[,"DIFSE"]      <- sqrt(TM[,"IFVar2"] + TM[,"IFVar1"])
  
    TS <- merge(TSA1,TSA2,by=c("Sample"))
    TS[,"Difference"] <- TS[,"IFEst2" ] - TS[,"IFEst1" ]
    TS[,"DIFSE"]      <- sqrt(TS[,"IFVar2"] + TS[,"IFVar1"])
    TS[,"TIF"]        <- (TS[,"Difference"] - TS[,"MIFEst2"] + TS[,"MIFEst1"]) / TS[,"DIFSE"]
  
    TS    <- data.frame(TS)

    ## quantiles
    cicut  <- CIlevel
    cicut1 <- (1 - CIlevel)/2
    cicut2 <- 1 - cicut1
    myq <- function(x) {
      x1 <- as.vector(quantile(x, c(cicut1,cicut2)))
      x2 <- as.vector(quantile(abs(x), c(cicut)))
      return(c(x1,x2))                
    }
##    ag1 <- matrix(myq( TS[,c("TVals")]), nrow=1 )
##    ag2 <- matrix(myq( TS[,c("Difference")]), nrow=1)
    ag3 <- matrix(myq( TS[,c("TIF")  ]), nrow=1)

##    colnames(ag1) <- c("tLow",   "tHigh",     "tSym")
##    colnames(ag2) <- c("DIFLow", "DIFHigh",   "DIFSym")
    colnames(ag3) <- c("tIFLow", "tIFHigh",   "tIFSym")

##    TM <- merge(TM,ag1)
##    TM <- merge(TM,ag2)
    TM <- merge(TM,ag3)

    zcut <- qnorm( cicut2 )
    ## IF CI
    lb1 <- TM[,"Difference"] - zcut * TM[,"DIFSE"]
    ub1 <- TM[,"Difference"] + zcut * TM[,"DIFSE"]

    ## Bootstrap CI
##    lb2 <- TM[,"Difference"] - 1.96 * TM[,"DIFSE" ]
##    ub2 <- TM[,"Difference"] + 1.96 * TM[,"DIFSE" ]

    ## jk CI
##    lb3 <- TM[,"Difference"] - 1.96 * TM[,"DjkSE" ]
##    ub3 <- TM[,"Difference"] + 1.96 * TM[,"DjkSE" ]

    ## quantiles
##    lb4 <- TM[,"DIFLow" ]
##    ub4 <- TM[,"DIFHigh"]

    ## t using IF se
    lb5 <- TM[,"Difference"] - TM[,"tIFHigh"] * TM[,"DIFSE"]
    ub5 <- TM[,"Difference"] - TM[,"tIFLow" ] * TM[,"DIFSE"]

    ## t using jk se
##    lb6 <- TM[,"Difference"] - TM[,"tHigh"] * TM[,"DjkSE"]
##    ub6 <- TM[,"Difference"] - TM[,"tLow" ] * TM[,"DjkSE"]

    ## t using IF se symmetric
    lb7 <- TM[,"Difference"] - TM[,"tIFSym"] * TM[,"DIFSE"]
    ub7 <- TM[,"Difference"] + TM[,"tIFSym"] * TM[,"DIFSE"]
  
    ## t using jk se symmetric
##    lb8 <- TM[,"Difference"] - TM[,"tIFSym"] * TM[,"DjkSE"]
##    ub8 <- TM[,"Difference"] + TM[,"tIFSym"] * TM[,"DjkSE"]

    Difference <- TM[,"Difference"]
    CI <- cbind( alpha1, alpha2,Difference,lb1, ub1, lb5, ub5, lb7, ub7)

    if ( i1 == 1 & i2 == 1) {
       CIout <- CI
       TMout <- TM
    }
    else {
       CIout <- rbind(CIout,CI)
       TMout <- rbind(TMout,TM)
    }
}
}
 Ret <- list( TM = TMout, CI = CIout, Nalpha = Nalpha, alphaList = alphaList, CIlevel=CIlevel )
 return(Ret)
}
