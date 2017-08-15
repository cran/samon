\name{samonDataCheck}
\alias{samonDataCheck}
\title{Samon Data Check}
\description{Checks data for intermittent missing or monotone missing pattern}
\usage{samonDataCheck( data )}
\arguments{
  \item{data}{ matrix with (i,j) entry representing value for subject i at
    time-point j. }
}

\details{
 Returns a basic summary of samon input data. This includes a check that
 there is not missing data at baseline. 
}
\value{
 Returns a list with items:  N the number or rows in data, NT the number
 of columns in data, missingBaseline which is set to 1 if any baseline
 data is missing, and, interMissing, indicating if there are any
 intermittent data.

 desc a matrix of dimension N by 4 with each row representing the same
 row in the original data.  Columns are:
 col 1: 1 if non-missing baseline value
 col 2: last available time
 col 3: last available observation
 col 4: number of non missing values

 Note that col2 - col4 gives the number of intermittent missing values.

 missingPatterns, a list of missing patterns found in the data
 NmissingTables, frequencies of missing patterns, and PmissingTable
 the proportion of each missing pattern in the data.
}
\examples{
 # monotone missing
 data("samonPANSS1")
 chkm <- samonDataCheck( samonPANSS1 )

 # non monotone missing
 data("VAS1")
 chknm <- samonDataCheck( VAS1 )
}
