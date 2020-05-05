\name{samonECompleterStatus}
\alias{samonECompleterStatus}
\title{Computes the difference in the expected value of non-completers and completers}
\description{Computes the difference in the expected value of
  non-completers and completers}
\usage{samonECompleterStatus( Y, BC )}
\arguments{
  \item{Y}{ vector of outcomes}
  \item{BC}{ a vector of estimates for the mean of Y}
}

\details{
 Computes the expected value for completers and non-completers and hence
 the difference between the two.
}
\value{
 Returns a vector of expected differences.
}
\examples{
data("samonPANSS1")
data( "PSummary1" )
PSummary1 <- (PSummary1$TM)[,c("alpha","IFEst")]
Pdiff1 <- samonECompleterStatus( samonPANSS1[,ncol(samonPANSS1)],
                                 PSummary1[,2] )

data("VAS1")
data("VSummary1" )
VSummary1 <- (VSummary1$TM)[,c("alpha","IFEst")]
Vdiff1 <- samonECompleterStatus( VAS1[,ncol(VAS1)], VSummary1[,2] )
}
