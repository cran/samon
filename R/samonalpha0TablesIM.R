alpha0TablesIM <- function( alpha0TableResultsObj, trt1lab = "Treatment 1", trt2lab = "Treatment 2" ) {

    Obsmeans <- alpha0TableResultsObj$Obsmeans

    ltrt1lab <- nchar(trt1lab)
    ltrt2lab <- nchar(trt2lab)

    lastp    <- round( (15 + 39 + ltrt1lab) / 2, 1)
    lastp2   <- round( ( 9 + 39 + ltrt2lab) / 2, 1)
    
    fmt      <- sprintf("%%%2.0fs %%%2.0fs",lastp,lastp2)
    
    dash    <- paste(rep("-",72), sep="", collapse="")
    headerT <- "                      Treatment 1                  Treatment 2        "
    headerT <- sprintf(fmt,trt1lab,trt2lab)
    header  <- "    Time        N       Mean         SD        N       Mean         SD"

    Text1 <- "Number of values available at each time-point, mean and standard"
    Text2 <- "deviation of outcome."
    cat("\n\n")
    cat(paste(Text1,  "\n"))
    cat(paste(Text2,  "\n"))
    cat(paste(dash,   "\n"))
    cat(paste(headerT,"\n"))
    cat(paste(header, "\n"))
    cat(paste(dash,   "\n"))

    for ( i in 1:nrow(Obsmeans) )  {
       line <-  sprintf("%8.0f %8.0f %10.2f %10.2f %8.0f %10.2f %10.2f\n", Obsmeans[i,1], Obsmeans[i,2], Obsmeans[i,3], Obsmeans[i,4], Obsmeans[i,5], Obsmeans[i,6], Obsmeans[i,7])
       cat(line)
    }    
    cat(paste(dash,  "\n\n\n\n"))


    # Table 2 obs
    ObsResults <- alpha0TableResultsObj$ObsResults
    ObsTtest   <- alpha0TableResultsObj$ObsTtest

    dash    <- paste(rep("-",75), sep="", collapse="")
    header  <- "Treatment                 N       Mean         SD    t-statistic    p-value"
    Text1 <- "At final time-point (observed)"
    cat(paste(Text1,  "\n"))
    cat(paste(dash,   "\n"))
    cat(paste(header, "\n"))
    cat(paste(dash,   "\n"))
    
    line    <- sprintf("%-16s %10.0f %10.2f %10.2f\n",trt1lab, ObsResults[1,"NComplete"],ObsResults[1,"meanComplete"],ObsResults[1,"sdComplete"])
    cat(line) 
    line    <- sprintf("%-16s %10.0f %10.2f %10.2f\n",trt2lab, ObsResults[2,"NComplete"],ObsResults[2,"meanComplete"],ObsResults[2,"sdComplete"])
    cat(line)

    cat("\n")
    line    <- sprintf("Difference (2-1)  %20.2f  %24.3f %10.3f\n", ObsTtest[1,"Difference"],ObsTtest[1,"ttestStatistic"],ObsTtest[1,"Pvalue"])
    cat(line)
    
    cat(paste(dash,   "\n\n\n"))



    # Table 3 LOCF
    LOCFResults <- alpha0TableResultsObj$LOCF
    LOCFTtest   <- alpha0TableResultsObj$LOCFTtest

    dash    <- paste(rep("-",75), sep="", collapse="")
    header  <- "Treatment                 N       Mean         SD    t-statistic    p-value"
    Text1 <- "At final time-point (LOCF)"
    cat(paste(Text1,  "\n"))
    cat(paste(dash,   "\n"))
    cat(paste(header, "\n"))
    cat(paste(dash,   "\n"))
    
    line    <- sprintf("%-16s %10.0f %10.2f %10.2f\n",trt1lab, LOCFResults[1,"N"],LOCFResults[1,"mean"],LOCFResults[1,"sd"])
    cat(line) 
    line    <- sprintf("%-16s %10.0f %10.2f %10.2f\n",trt2lab, LOCFResults[2,"N"],LOCFResults[2,"mean"],LOCFResults[2,"sd"])
    cat(line)

    cat("\n")
    line    <- sprintf("Difference (2-1)  %20.2f  %24.3f %10.3f\n", LOCFTtest[1,"Difference"],LOCFTtest[1,"ttestStatistic"],LOCFTtest[1,"Pvalue"])
    cat(line)
    
    cat(paste(dash,   "\n\n\n"))


    # Table 4 IF results
    IFResults <- alpha0TableResultsObj$IFResults

    CIlevel  <- 100*alpha0TableResultsObj$CIlevel
    if ( round(CIlevel) == CIlevel ) {
      headerT1 <- sprintf("                           IF            Bootstrap %2d%%           Jackknife %2d%%", CIlevel, CIlevel)
      headerT1 <- sprintf("                                IF              MI Bootstrap %2d%%", CIlevel )
    } else {
      headerT1 <- sprintf("                           IF          Bootstrap %4.1f%%         Jackknife %4.1f%%", CIlevel, CIlevel)
      headerT1 <- sprintf("                                IF            MI Bootstrap %4.1f%%", CIlevel )
    }
        
    dash    <- paste(rep("-",64), sep="", collapse="")
    headerT <- "                                IF                 MI Bootstrap "
    header  <- "Treatment                    Estimate             LCL        UCL"
    Text1 <- "At final time-point alpha = 0."
    cat(paste(Text1,  "\n"))
    cat(paste(dash,   "\n"))
##    cat(paste(headerT,"\n"))
    cat(paste(headerT1,"\n"))
    cat(paste(header, "\n"))
    cat(paste(dash,   "\n"))
    
    line    <- sprintf("%-16s %20.3f %15.3f %10.3f\n",trt1lab, IFResults[1,"Estimate"],IFResults[1,"LCLBootstrap"],IFResults[1,"UCLBootstrap"])
    cat(line) 
    line    <- sprintf("%-16s %20.3f %15.3f %10.3f\n",trt2lab, IFResults[2,"Estimate"],IFResults[2,"LCLBootstrap"],IFResults[2,"UCLBootstrap"])
    cat(line)

    cat("\n")
    line    <- sprintf("Difference (2-1)  %19.3f %15.3f %10.3f\n", IFResults[3,"Estimate"],IFResults[3,"LCLBootstrap"],IFResults[3,"UCLBootstrap"])
    cat(line)
    
    cat(paste(dash,   "\n\n\n"))



    # Table 5 Optimal sigmaH and sigmaF
    optRes <- alpha0TableResultsObj$optRes

    dash    <- paste(rep("-",39), sep="", collapse="")
    header  <- "Treatment            SigmaH     SigmaF"
    Text1 <- "Smoothing parameters:"
    cat(paste(Text1,  "\n"))
    cat(paste(dash,   "\n"))
    cat(paste(header, "\n"))
    cat(paste(dash,   "\n"))
    
##    line    <- sprintf("%-16s %10.3f %10.3f\n",trt1lab, optRes[1,"SigmaP"],optRes[1,"SigmaQ"])
##    cat(line) 
##    line    <- sprintf("%-16s %10.3f %10.3f\n",trt2lab, optRes[2,"SigmaP"],optRes[2,"SigmaQ"])
##    cat(line)

    for ( i in 1:nrow(optRes) )  {
       if ( i == 1 ) nm <- trt1lab
       else if ( i == ((nrow(optRes)/2)+1)) nm <- trt2lab
       else nm <- " "         
       line    <- sprintf("%-16s %10.3f %10.3f\n",nm, optRes[i,"SigmaH"],optRes[i,"SigmaF"])
       cat(line)
    }    
##    cat(paste(dash,  "\n\n\n\n"))
    cat(paste(dash,   "\n\n"))

}
