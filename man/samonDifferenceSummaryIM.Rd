\name{samonDifferenceSummaryIM}
\alias{samonDifferenceSummaryIM}
\title{Summary function for the difference in treatment}
\description{Produces a summary of the difference in treatment effect
  using two samonIM objects. It combines bootstrap imputed results
  to produce confidence intervals}
\usage{samonDifferenceSummaryIM( trt1, trt2, CIlevel = 0.95 )}
\arguments{
  \item{trt1}{ the result from a call to samonSummaryIM }
  \item{trt2}{ the result from a call to samonSummaryIM }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines estimates (imputed) bootstrap estimates to  produce summary
 estimates and confidence intervals for an estimate of the difference in
 treatment effect from two samon objects.

 The difference in treatment estimates is taken to be those of treatment
 2 minus those of treatment 1 at the same value of sensitivity parameter
 alpha.
}
\value{
 Returns a list. Items include: TM which contains treatment estimates
 and the difference, together with standard error estimates based on
 multiple imputation of bootstrap samples. It also contains CI which
 contains confidence intervals. Other items are returned for convenience
 are n10, the number of subjects in the first treatment group, and n20
 the number of subjects in the second treatment group.
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
