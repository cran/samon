\name{samonCombineIM}
\alias{samonCombineIM}
\title{Combines results from separate runs of the samonIM function}
\description{Combines lists returned from the samonIM function to form a
  single object.}
\usage{samonCombineIM( filenames, replaceSampleNo=TRUE, ForceAppend =
FALSE, objlist = NULL )}
\arguments{
  \item{filenames}{ a list of rds filenames each containing the output
    from samonIM}
  \item{replaceSampleNo}{ A logical value determining if the sample
    numbers should be replaced in the output object.  Usually the sample
  numbers are 1,2,...,NSamples in each file. When the results in
  filenames are from different runs of samon with different seeds, the
  results in the output are given unique sample numbers, usually 1, 2,
  ... NTot, where NTot is the total number of samples.}
  \item{ForceAppend}{ put results together even if the files don't match
    on some criteria, such as the start values for sigma H or sigma
    F. Otherwise report the problem and return.}
  \item{objlist}{ an alternative input of samon objects to be combined.}
}
\details{
 When the samonIM function is used in parallel a number of results
 are returned.  If these results are stored in rds files, the
 samonCombine function can extract the results and combine them to
 produce a single object.
 
 This is usually done for a single treatment arm at a time.
}
\value{
 Returns a list similar to that returned by samonIM itself.  This list may
 be passed to samonSummaryIM and other samon functions. 
}
\examples{
# outputs from samonIM -- run on VAS1 with different seeds
data("V1Res001")
data("V1Res002")

# combine them
V1Results <- samonCombineIM( objlist = list(V1Res001,V1Res002) )
}
