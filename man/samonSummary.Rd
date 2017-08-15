\name{samonSummary}
\alias{samonSummary}
\title{Summary function for an object returned by the samon function}
\description{Produces summaries of samon objects, combining bootstrap and
  jackknife results to produce confidence intervals}
\usage{samonSummary( trt, CIlevel = 0.95 )}
\arguments{
  \item{trt}{ the result from a call to samon or samonCombine }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines estimates, bootstrap estimates and jackknife estimates to
 produce summary estimates and confidence intervals for a samon object.
}
\value{
 Returns a list. Items include TM which contains treatment estimates and
 standard errors for the main data; TS contains estimates for the
 bootstrap samples and includes jackknifes standard error estimates;
 CI contains confidence intervals. 
}
\examples{
data("P1Results")
Summary1 <- samonSummary( P1Results )
}
