\name{samonCrossSummaryIM}
\alias{samonCrossSummaryIM}
\title{Summary function for the difference in effect of two treatments
  at different sensitivity parameters.}
\description{Produces a summary of the difference in treatment effect
  using two samonSummaryIM objects, Differences and confidence intervals
  are produced for each pair of sensitivity parameters.}
\usage{samonCrossSummaryIM( trt1, trt2, CIlevel = 0.95 )}
\arguments{
  \item{trt1}{ the result from a call to samonSummaryIM for treatment 1 }
  \item{trt2}{ the result from a call to samonSummaryIM for treatment 2 }
  \item{CIlevel}{ the confidence level for confidence intervals }
}

\details{
 Combines bootstrap estimates to produce summary estimates and
 confidence intervals for the difference in  estimates from two
 treatment groups.  The difference is computed for all combinations of
 the smoothing parameter, alpha, found in the Sum1 input object. 

 The difference in treatment estimates is taken to be those of treatment
 2 minus those of treatment 1.
}
\value{
 Returns a list. This list includes TM which contains treatment estimates,
 their difference, and, standard error estimates. The list also contains,
 CI, which contains confidence intervals.
 Other items are returned for convenience are n10, the number of
 subjects in the first treatment group, and n20 the number of subjects
 in the second treatment group.
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
