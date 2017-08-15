\name{samonDifferenceSummary}
\alias{samonDifferenceSummary}
\title{Summary function for the difference in treatment effect}
\description{Produces a summary of the difference in treatment effect
  using two samon objects. It combines bootstrap and jackknife results
  to produce confidence intervals}
\usage{samonDifferenceSummary( trt1, trt2, CIlevel = 0.95 )}
\arguments{
  \item{trt1}{ the result from a call to samonSummary }
  \item{trt2}{ the result from a call to samonSummary }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines estimates, bootstrap estimates and jackknife standard error
 estimates to produce summary estimates and confidence intervals for an
 estimate of the difference in treatment effect from two samon objects.

 The difference in treatment estimates is taken to be those of treatment
 2 minus those of treatment 1 at the same value of sensitivity parameter
 alpha.
}
\value{
 Returns a list. Items include: TM which contains treatment estimates,
 standard errors for these estimates, the estimate of their difference
 and estimates of the standard error of the difference.  Also contains
 CI which contains  confidence intervals.
 Other items are returned for convenience are n10, the number of
 subjects in the first treatment group, and n20 the number of subjects
 in the second treatment group.
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
