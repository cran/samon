\name{alpha0Tables}
\alias{alpha0Tables}
\alias{alpha0TableResults}
\title{Makes samon result tables at alpha = 0}
\description{alpha0TableResults collects results for tables at sensitivity parameter
  alpha = 0 and alpha0Tables prints the output from alpha0TableResults} 
\usage{alpha0TableResults( Y1, samonR1, Y2, samonR2, samonRD)
alpha0Tables( alpha0TableResultsObj, trt1lab = "Treatment 1",
trt2lab = "Treatment 2")}
\arguments{
  \item{Y1}{ treatment 1's raw data.}
  \item{samonR1}{ samon results for treatment 1.}
  \item{Y2}{ treatment 2's raw data.}
  \item{samonR2}{ samon results for treatment 2.}
  \item{samonRD}{ samon results for the difference of treatment 2 minus
    treatment 1.}
  \item{alpha0TableResultsObj}{ a list produced by the function alpha0TableResults.}
  \item{trt1lab}{ a label for treatment 1.}
  \item{trt2lab}{ a label for treatment 2.}
}

\details{
 alpha0TableResults computes material to build the table and
 alpha0Tables prints tables as simple text.
}
\value{
 alpha0TableResults returns table data. 
 alpha0Tables returns nothing.
}
\examples{
data("samonPANSS1")
data("samonPANSS2")

data("PSummary1")
data("PSummary2")
data("PSummaryD")

alpha0Results <- alpha0TableResults( samonPANSS1, PSummary1,
                                     samonPANSS2, PSummary2,
                                                  PSummaryD )
alpha0Tables(alpha0Results)
}

