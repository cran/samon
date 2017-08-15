\name{The VAS trial.}
\docType{data}
\alias{VAS1}
\alias{VAS2}
\title{The VAS randomized trial.}
\description{These dataframes store data from a repeated measures
  clinical trial.  VAS1 contains data from the control arm and VAS2
  contains data from the active arm.}
\details{In each of VAS1 and VAS2 a row represtents the outcome measure
  (VAS) for one subject at each of 9 timepoints. 
  Both treatment arms contain intermittent missing values.
  \itemize{
  \item V1 VAS score at baseline 
  \item V2 VAS score at follow-up time-point 1
  \item V3 VAS score at follow-up time-point 2
  \item ...
  \item V9 VAS score at follow-up time-point 8 the last time-point.
}
}
\usage{
data(VAS1)
data(VAS2)
}
\format{
  VAS1 is a 255 by 9 dataframe
  VAS2 is a 256 by 9 dataframe
}

