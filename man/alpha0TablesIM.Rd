\name{alpha0TablesIM}
\alias{alpha0TablesIM}
\alias{alpha0TableResultsIM}
\title{Makes samonIM result tables at sensitivity parameter alpha = 0.}
\description{alpha0TableResultsIM collects results for tables at
  sensitivity parameter alpha = 0 and alpha0TablesIM prints the output
  from alpha0TableResultsIM} 
\usage{alpha0TableResultsIM( Y1, samonR1, Y2, samonR2, samonRD)
alpha0TablesIM( alpha0TableResultsObj, trt1lab = "Treatment 1",
trt2lab = "Treatment 2")}
\arguments{
  \item{Y1}{ treatment 1's raw data.}
  \item{samonR1}{ samon results for treatment 1.}
  \item{Y2}{ treatment 2's raw data.}
  \item{samonR2}{ samon results for treatment 2.}
  \item{samonRD}{ samon results for the difference of treatment 2 minus
    treatment 1.}
  \item{alpha0TableResultsObj}{ a list produced by the function
    alpha0TableResultsIM.}
  \item{trt1lab}{ a label for treatment 1.}
  \item{trt2lab}{ a label for treatment 2.}
}
\details{
 alpha0TableResultsIM computes material to build the table and
 alpha0TablesIM prints tables as simple text.
}
\value{
 alpha0TableResultsIM returns table data. 
 alpha0TablesIM returns nothing.
}
\examples{
data("VAS1")
data("VAS2")

data("VSummary1")
data("VSummary2")
data("VSummaryD")

alpha0Results <- alpha0TableResultsIM( VAS1, VSummary1,
                                       VAS2, VSummary2,
                                             VSummaryD )
alpha0TablesIM(alpha0Results)
}
