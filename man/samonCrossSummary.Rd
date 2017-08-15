\name{samonCrossSummary}
\alias{samonCrossSummary}
\title{Summary function for the difference in effect of two treatments
  at different sensitivity parameters alpha.}
\description{Produces a summary of the difference in treatment effect
  using two samonSummary objects. Differences and confidence intervals
  are produced for each pair of sensitivity parameters.}
\usage{samonCrossSummary( trt1, trt2, CIlevel = 0.95 )}
\arguments{
  \item{trt1}{ the result from a call to samonSummary }
  \item{trt2}{ the result from a call to samonSummary }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines estimates, bootstrap estimates and jackknife estimates to
 produce summary estimates and confidence intervals for the difference
 in  estimates between two treatment groups.  The difference is computed
 for all combinations of the smoothing parameter alpha found in the Sum1
 input object. 

 The difference in treatment estimates is taken to be those of treatment
 2 minus those of treatment 1.
}
\value{
 Returns a list. Items include TM which contains treatment estimates,
 their difference, and, standard error estimates and CI which contains
 confidence intervals. Other items are returned for convenience are n10,
 the number of subjects in the first treatment group, and n20 the number
 of subjects in the second treatment group.
}
\examples{
# samon results for treatment 1 and 2 respectively:
data("P1Results")
data("P2Results")

Summary1 <- samonSummary( P1Results )
Summary2 <- samonSummary( P2Results )
SummaryD <- samonDifferenceSummary( Summary1, Summary2 )
SummaryX <- samonCrossSummary( Summary1, Summary2 )
}
