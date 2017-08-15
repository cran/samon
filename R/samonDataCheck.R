# Samon: Summarizes data in an input data matrix
# particularly summarizes missing data patterns.
# ----------------------------------------------------
samonDataCheck <- function( data ) {
   NT <- ncol(data)
   N  <- nrow(data)

   min = min(data, na.rm=TRUE)
   max = max(data, na.rm=TRUE)

   ## desc[,1] Valid (non-missing) baseline 
   ## desc[,2] last available t
   ## desc[,3] last available value
   ## desc[,4] number of observed values
   desc     <- matrix(0,N,4)
   missPattern <- rep(paste( rep("*",NT), sep="",collapse=""),N)
   for ( i in 1:N  ) {
       desc[i,]  <- rep(0,4)
       for ( t in 1:NT ) {
           if ( is.na(data[i,t]))  substr(missPattern[i],t,t+1) <- "_"
           if ( (!is.na(data[i,t]))) {
               desc[i,2] <- t
               desc[i,3] <- data[i,t]
               desc[i,4] <- desc[i,4] + 1
           }
       }
   }
   desc[,1] <- 1 - is.na(data[,1])

   missingBaseline     <- any(desc[,1] == 0)
   NinterMiss          <- desc[,2] - desc[,4] 
   intermittentMissing <- any( NinterMiss != 0 )
   completeData        <- sum( desc[,4] == NT )

   cat("\n\n")
   cat("Samon Data Check:\n")
   cat("--------------------------------------------------\n")
   cat(sprintf("Number of timepoints:                   %9.0f\n",                NT  ))
   cat(sprintf("Number of subjects:                     %9.0f\n",                 N  ))
   cat(sprintf("Minimum observed value:                 %9.0f\n",               min  ))
   cat(sprintf("Maximum observed value:                 %9.0f\n",               max  ))
   cat(sprintf("Average number of timepoints on study:  %9.2f\n",    mean(desc[,2])  ))
   cat(sprintf("Total number of observed values:        %9.0f\n",     sum(desc[,4])  ))
   cat(sprintf("Subjects observed at final timepoint:   %9.0f\n", sum(desc[,2]==NT)  ))
   cat(sprintf("Subjects observed at all timepoints:    %9.0f\n",      completeData  ))
   cat("\n")

   if ( missingBaseline == 1 ) {    
       cat(sprintf("Missing baseline data found:\n"))
       cat(sprintf("   subjects = %8.0f \n", N - sum(data[,1])))
       cat("\n")
   }
   if ( intermittentMissing == 1 ) {     
       cat(sprintf("Intermittent Missing data found:\n"))
       cat(sprintf("    subjects with IM    = %9.0f \n", N - sum( desc[,2] == desc[,4] )))    
       cat(sprintf("    number of IM values = %9.0f \n", sum( desc[,2] - desc[,4] )))    
       cat("\n")
   }

   ntab <- table(missPattern)
   ptab <- prop.table(ntab)

   missHead = paste( paste( rep(" ",NT+3), sep="", collapse=""), paste( rep(" ",7), sep="", collapse=""), "N", paste( rep(" ",1), sep="", collapse=""), "proportion", "\n")
   cat("\n")
   cat("Missing Patterns:\n")
   cat(missHead)
   for ( nm in names(ntab) ) {
       cat( nm, " : ", format(ntab[nm], justify="right", width=8), format(round(ptab[nm],4), nsmall=4, digits=4, justify="right", width=12), "\n" )
   }
   cat("\n\n")

   dimnames(desc)[[2]] <- list("baseline","lastTime","lastValue","NValues")
   ret <- list( N = N, NT = NT, missingBaseline = missingBaseline, intermittentMissing = intermittentMissing, desc = desc, missingPatterns = missPattern, NmissingTable=ntab, PmissingTable=ptab )
   return(ret)
}
