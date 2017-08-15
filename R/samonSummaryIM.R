# --------------------------------------------------------------------------------
# SummaryIM function for samon objects
# --------------------------------------------------------------------------------
samonSummaryIM <- function(trt, CIlevel=0.95) {

  trtnames  <- names(trt)

  NFills    <- trt$Nfills
  M         <- trt$NIMimpute
 
  HM        <- trt$HM   [ trt$HM [ ,2] <= M, c("Type","Convergence","Iterations","SigmaH","lossH")]
  FM        <- trt$FM   [ trt$FM [ ,2] <= M, c("Type","Convergence","Iterations","SigmaF","lossF")]
  TM        <- trt$IFM  [ trt$IFM[ ,2] <= M, c("Type","alpha","AEst","AVar","IFEst","IFVar")]  
  TS        <- trt$IFS  [ trt$IFS[ ,2] <= M, c("Sample","Type","alpha","AEst","AVar","IFEst","IFVar")]
  
  Nalpha    <- trt$Nalpha
  alphaList <- trt$alphaList
  NSamples  <- trt$NSamples
  n0        <- trt$n0

  HM    <- data.frame(HM)
  HMagg <- aggregate( . ~ 1, HM[,c("Convergence","Iterations","SigmaH","lossH")], function(x) { c(mean(x),var(x),min(x),max(x)) })
  HMagg <- matrix(unlist(HMagg),4,byrow=TRUE,dimnames=list(c("Convergence","Iterations","SigmaH","lossH"),c("mean","sd","min","max")))
  
  FM    <- data.frame(FM)
  FMagg <- aggregate( . ~ 1, FM[,c("Convergence","Iterations","SigmaF","lossF")], function(x) { c(mean(x),var(x),min(x),max(x)) })
  FMagg <- matrix(unlist(FMagg),4,byrow=TRUE,dimnames=list(c("Convergence","Iterations","SigmaF","lossF"),c("mean","sd","min","max")))

  TM    <- data.frame(TM)
  TMagg <- aggregate( . ~ alpha, TM[,c("alpha","AEst","IFEst","IFVar")], function(x) { c(mean(x),var(x)) } )
  
  TMagg <- cbind(TMagg[[1]],TMagg[[2]],TMagg[[3]],TMagg[[4]])
  TMagg <- data.frame(TMagg)
  colnames(TMagg) <- c("alpha","mAEst","vAEst","mIFEst","vIFEst","mIFVar","vIFVar")
  TMagg[,"MIIFVar"] <- TMagg[,"mIFVar"] + ( 1 + 1/M ) * TMagg[,"vIFEst"] 
  TMagg <- TMagg[,c("alpha","mAEst","mIFEst","MIIFVar")]
  colnames(TMagg) <- c("alpha","AEst","IFEst","MIIFVar")

  TS    <- data.frame(TS)
  TSagg <- aggregate( . ~ Sample + alpha, TS[,c("Sample","alpha","IFEst","IFVar")], function(x) { c(mean(x),var(x)) } )
  TSagg <- cbind(TSagg[[1]],TSagg[[2]],TSagg[[3]],TSagg[[4]])
  TSagg <- data.frame(TSagg)
  colnames(TSagg) <- c("Sample","alpha","mIFEst","vIFEst","mIFVar","vIFVar")
  TSagg[,"MIIFVar"] <- TSagg[,"mIFVar"] + ( 1 + 1/M ) * TSagg[,"vIFEst"]
  TSagg <- TSagg[,c("Sample","alpha","mIFEst","MIIFVar")]
  colnames(TSagg) <- c("Sample","alpha","IFEst","MIIFVar")

  ## add main means to bootstraps
  TMSelect <- TMagg[,c("alpha","IFEst")]
  colnames(TMSelect) <- c("alpha","mIFEst")

  TSagg <- merge(TSagg,TMSelect,by=c("alpha"))

  ## calculate a couple of t values
  TSagg[,"TIF"  ] <- ( TSagg[,"IFEst"] - TSagg[,"mIFEst"] ) / sqrt(TSagg[,"MIIFVar"])

  ## quantiles
  cicut  <- CIlevel
  cicut1 <- (1 - CIlevel)/2
  cicut2 <- 1 - cicut1
  myq <- function(x) {
    x1 <- as.vector(quantile(x, c(cicut1,cicut2)))
    x2 <- as.vector(quantile(abs(x), c(cicut)))
    return(c(x1,x2))                
  }
  ag3 <- aggregate( . ~ alpha, TSagg[,c("alpha","TIF")   ], myq )
  ag3 <- cbind(ag3[1],ag3[[2]])
  
  colnames(ag3) <- c("alpha", "tIFLow", "tIFHigh", "tIFSym")

  TM <- merge(TMagg,ag3,by=c("alpha"))

  zcut <- qnorm( cicut2 )
  ## IF CI
  lb1 <- TM[,"IFEst"] - zcut * sqrt(TM[,"MIIFVar"])
  ub1 <- TM[,"IFEst"] + zcut * sqrt(TM[,"MIIFVar"])

  ## t using IF SE asymmetric
  lb5 <- TM[,"IFEst"] - TM[,"tIFHigh"] * sqrt(TM[,"MIIFVar"])
  ub5 <- TM[,"IFEst"] - TM[,"tIFLow" ] * sqrt(TM[,"MIIFVar"])

  ## t using IF se symmetric
  lb7 <- TM[,"IFEst"] - TM[,"tIFSym"] * sqrt(TM[,"MIIFVar"])
  ub7 <- TM[,"IFEst"] + TM[,"tIFSym"] * sqrt(TM[,"MIIFVar"])

  CI <- cbind(alphaList, lb1, ub1, lb5, ub5, lb7, ub7)
  Ret <- list( TM = TM, TS = TSagg, CI = CI, n0 = n0, NSamples = NSamples, Nalpha = Nalpha, alphaList = alphaList, HM = HM, HMSummary = HMagg, FM = FM, FMSummary = FMagg, CIlevel=CIlevel )
  return(Ret)
}
