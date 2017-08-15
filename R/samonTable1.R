## -----------------------------------------------------
## Display a matrix produced by the Tabmat1 function
## -----------------------------------------------------
samonTable1 <- function( Tabdata, trtlab = "Treatment 1" ) {
    
  dash     <- paste(rep("-",108), sep="", collapse="")

  headerT2 <- "                                     Last Observed               Intermittent Missing"
  headerT3 <- "                            ------------------------------     -----------------------"
  headerT4 <- "                                      As % of     As % of                     As % of";
  headerT5 <- "          On-         N                Number      Number                      Number            Observed    "
  headerT6 <- "   T    Study       Obs         N    On-Study    Observed            N       On-Study         mean        SD "
  
  Text1 <- trtlab
  Text2 <- "Number of observed values at each time-point, mean and standard deviation of observed outcome."
  cat("\n\n")
  cat(paste(Text1,   "\n"))
  cat(paste(Text2,   "\n"))
  cat(paste(dash,    "\n"))
  cat(paste(headerT2,"\n"))
  cat(paste(headerT3,"\n"))
  cat(paste(headerT4,"\n"))
  cat(paste(headerT5,"\n"))
  cat(paste(headerT6,"\n"))
  cat(paste(dash,    "\n"))
  
  for ( i in 1:nrow(Tabdata) )  {
      if ( i == 1 ) {
          line <-  sprintf("%4.0f %8.0f %9.0f %9.0f %11.2f %11.2f                            ", Tabdata[i,1], Tabdata[i,2], Tabdata[i,3], Tabdata[i,4], 100*Tabdata[i,5], 100*Tabdata[i,6])
      } else if ( i == nrow(Tabdata) ) {
          line <-  sprintf("%4.0f %8.0f %9.0f %9.0f %11.2f %11.2f                            ", Tabdata[i,1], Tabdata[i,2], Tabdata[i,3], Tabdata[i,4], 100*Tabdata[i,5], 100*Tabdata[i,6])
      } else {
          line <-  sprintf("%4.0f %8.0f %9.0f %9.0f %11.2f %11.2f %12.0f %14.2f", Tabdata[i,1], Tabdata[i,2], Tabdata[i,3], Tabdata[i,4], 100*Tabdata[i,5], 100*Tabdata[i,6], Tabdata[i,7], 100*Tabdata[i,8])
      }
      cat(line)
      if ( i == 1 || i == nrow(Tabdata) ) Ndo = 1
      else Ndo = 1
      for ( j in 1:Ndo ) {
          line <-  sprintf("%13.3f %9.3f", Tabdata[i,8 + (j-1) * 2 + 1], Tabdata[i,8 + (j-1) * 2 + 2])
          cat(line)
      }
      cat("\n")
  }    
  cat(paste(dash,  "\n\n\n\n"))
}    


