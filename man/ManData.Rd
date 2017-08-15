\name{ManData}
\alias{P1Res001}
\alias{P1Res002}
\alias{P1Results}
\alias{P2Results}
\alias{PSummary1}
\alias{PSummary2}
\alias{PSummaryD}
\alias{V1Res001}
\alias{V1Res002}
\alias{V1Results}
\alias{V2Results}
\alias{VSummary1}
\alias{VSummary2}
\alias{VSummaryD}
\docType{data}
\title{Data used in the man page examples.}
\description{Some functions within the 'samon' package take some time to
  run, particularly when the number of bootstraps requested is large.
  To illustrate some of the functions on the man pages pre-packaged
  results are sometime used to make the examples run faster. Here we
  list some of this data.
}
\details{The samonPANSS1 and samonPANSS2 datasets are documented on
  their own man pages.

 \describe{
  \item{P1Res001 and P1Res002}{ are two small output objects from the
  samon function and are used to illustrate the functioning of the
  samonCombine function.}
  \item{P1Results and P2Results}{ are combined results for samonPANSS1
  and samonPANSS2 respectively and are used to illustrate the
  functioning of the samonSummary function.}
  \item{PSummary1 and PSummary2}{ are summary objects from the
  samonSummary function and are used to illustrate the use of the
  samonDifferenceSummary and samonCrossSummary functions.}
  \item{PSummryD}{ this is the results from the samonDifferenceSummary
  function and is used along with PSummary1 and Psummary2 to illustrate
  the alpha0Table and samonEYCompleterStatus functions.}

  The VAS1 and VAS2 datasets are documented on their own man pages.

  \item{V1Res001 and V1Res002}{ are two small output objects returned
  from the samonIM function and are used to illustrate the functioning
  of the samonCombineIM function.}
  \item{V1Results and V2Results}{ are combined results for VAS1 and VAS2
  respectively and are used to illustrate the functioning of the
  samonSummaryIM function.}
  \item{VSummary1 and VSummary2}{ are summary objects from the
  samonSummaryIM function and are used to illustrate the use of the
  samonDifferenceSummaryIM and samonCrossSummaryIM functions.}
  \item{VSummryD}{ is the results from the samonDifferenceSummaryIM
  function and is used along with VSummary1 and VSummary2 to illustrate
  the alpha0Table and samonEYCompleterStatus functions.}
}}
\source{Pre-caluclated results from packaged data}
