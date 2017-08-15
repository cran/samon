\name{samonSummaryIM}
\alias{samonSummaryIM}
\title{Summary function for an object returned by the samonIM function}
\description{Produces summaries of samonIM objects, combining multiple
  imputed results to produce confidence intervals}
\usage{samonSummaryIM( trt, CIlevel = 0.95 )}
\arguments{
  \item{trt}{ the result from a call to samonIM or samonCombineIM }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines multiple imputation estimates and multiple imputation
 bootstrap estimates to produce summary estimates and confidence
 intervals for a samonIM object.
}
\value{
 Returns a list. Items include TM which contains treatment estimates and
 multiple imputation standard errors for the main data; TS contains
 bootstrap estimates and multiple imputation standard error bootstrap estimates.
 CI contains confidence intervals. 
}
\examples{
# V1Results and V2Results are the returned objects from the
# samonIM function associated with VAS treatment 1 and 2 respectively.
data("V1Results")
data("V2Results")

# summarize each arm, their difference and their cross difference.
VSummary1 <- samonSummaryIM(V1Results)
VSummary2 <- samonSummaryIM(V2Results)
VSummaryD <- samonDifferenceSummaryIM(VSummary1,VSummary2)
VSummaryX <- samonCrossSummaryIM(VSummary1,VSummary2)
}
