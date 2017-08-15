\name{The  PANSS trial.}
\docType{data}
\alias{samonPANSS1}
\alias{samonPANSS2}
\title{The PANSS randomized trial.}
\description{These dataframes store data from a repeated measures
  clinical trial.  samonPANSS1 contains data from the placebo arm and
  samonPANSS2 contains data from the active arm.}
\details{In each of samonPANSS1 and samonPANSS2 a row represents the
  outcome meansure (PANSS) for one subject at each of 6 timepoints.
  \itemize{
  \item V1 PANSS score at baseline 
  \item V2 PANSS score at follow-up time-point 1
  \item V3 PANSS score at follow-up time-point 2
  \item ...
  \item V6 PANSS score at follow-up time-point 5 the last time-point.
}
}
\usage{
data(samonPANSS1)
data(samonPANSS2)
}
\format{
  samonPANSS1 is a 88 by 6 dataframe
  samonPANSS2 is a 86 by 6 dataframe
}
\source{Data are simulated for the samon package}

