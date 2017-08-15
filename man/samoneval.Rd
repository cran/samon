\name{samoneval}
\alias{samoneval}
\title{Evaluates the loss function on a range of sigmas}
\description{Computes the loss function for a range of sigma}
\usage{samoneval( mat, Npart = 10, sigmaList = c(1), type="both")}
\arguments{
  \item{mat}{ matrix with (i,j) entry representing value for subject i at
    time-point j. }
  \item{Npart}{ Data is partitioned into Npart parts when estimating
    lossH and lossF.}
  \item{sigmaList}{ vector of sigmas at which to evaluate the loss function. }
  \item{type}{ one of "H", "F" or "both" to evaluate lossH, lossF or
    both for sigmaH or sigmaF in indicated range. }
}

\details{
 Stepping through the values in sigmaList this functions computes the
 associated loss function.  If type is set to "H" or "h", the loss
 function associated with dropout is computed.  If type is set to "F" or
 "f" the loss function associated with outcome is computed.  If type
 is set any other value then both types of loss function are evaluated.
}
\value{
 Returns an N by 2 matrix of results if type is "H" or "F" or an N by 4
 matrix if type is "both". There is one row for each value of sigma. If
 type is "H" then the columns are sigmaH and lossH, if type is "F" then
 the columns are sigmaF and lossF, and, if type = "both" then the
 columns are sigmaH, lossH, sigmaF, and lossF.
}
\examples{
data("samonPANSS1")
## dropout loss function
H1 <- samoneval( mat = samonPANSS1, Npart = 5,
      sigmaList      = seq(1,10,by=1),          
      type           = "H" )
}
