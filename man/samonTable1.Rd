\name{samonTable1}
\alias{samonTable1}
\alias{samonTabmat1}
\title{Summarizes two treatment groups.}
\description{samonTabmat1 takes a dataset Y and produces simple summary
  statistics which can be printed with the samontTable1 function.}
\usage{samonTabmat1( Y )
samonTable1( Tabdata, trtlab = "Treatment 1")}
\arguments{
  \item{Y}{ Data frame or matrix of data.  Rows are individual
    measuremnets and columns represent time-points. }
  \item{Tabdata}{ Matrix of summary stats produced by the samonTabmat1 function.}
  \item{trtlab}{ Label to place above table. }
}
\value{
samonTabmat1 returns an NT by 10 matrix where NT is the number of
columns in the input matrix Y.  The output matrix contains the columns
    1. t the time-point.
    2. Number On Study
    3. Number Observed
    4. Number last seen at time t
    5. Proportion last seen (of number on-study) 
    6. Proportion last seen (of number observed)
    7. Number intermittent missing data
    8. proportion intermittent missing (of number on-study)
    9. mean observed value at time t
   10. std observed value at time t

  samonTable1 does not return any object.  
}
\examples{
 # monotone missing
 data("samonPANSS1")
 samonTable1( samonTabmat1( samonPANSS1 ), trtlab = "PANSS Placebo" )

 # non monotone missing
 data("VAS1")
 samonTable1( samonTabmat1( VAS1 ), trtlab = "VAS Placebo" )
}


