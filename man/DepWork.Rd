\name{The DepWork trial}
\docType{data}
\alias{DepWork1}
\alias{DepWork2}
\title{The depression and working status trial.}
\description{These dataframes store data from a placebo controlled
  repeated measures clinical trial. DepWork1 contains data from the
  control arm and DepWork2 contains data from the active arm.}
\details{In each of DepWork1 and DepWork2 a row represents the outcome
  measure (working status) for one subject at each of 4 timepoints.
  At a given timepoint, working status is coded as 1 if an individual
  is not employed at that timepoint, and, as 2 if the individual is
  employed.  No individual was working at baseline.
  Both treatment arms contain intermittent missing values.
  \itemize{
  \item V1 work status at baseline 
  \item V2 work status at follow-up time-point 1
  \item V3 work status at follow-up time-point 2
  \item V4 work status at follow-up time-point 3
}
}
\usage{
data(DepWork1)
data(DepWork2)
}
\format{
  DepWork1 is a   552 by 4 dataframe
  DepWork2 is a  1249 by 4 dataframe
}

