\name{samonevalIM}
\alias{samonevalIM}
\title{Evaluates the loss function at a range of smoothing parameters}
\description{Computes the loss function for a range of smoothing parameters}
\usage{samonevalIM( mat, Npart = 10, sigmaList = c(1), inmodel =
inmodel, seed = 1, type="both")}
\arguments{
  \item{mat}{ matrix with (i,j) entry representing value for subject i at
    time-point j. }
  \item{Npart}{ number of parts to partition the data when evaluating
    the loss functions lossH and lossF. }
  \item{sigmaList}{ vector of sigmas at which to evaluate the loss function. }
  \item{inmodel}{ NT by 6 matrix indicating which variables to include
    in logistic model. }
  \item{seed}{ seed to use when imputing intermittent data. }
  \item{type}{ one of "H", "F" or "both" to evaluate lossH, lossF or
    both for sigmaH or sigmaF in indicated range. }
}

\details{
 Stepping through the values in sigmaList this functions computes the
 values of a loss function.  If type is set to "H" or "h" then the loss
 function associated with dropout is computed.  If type is set to "F" or
 "f" then the loss function associated with outcome is computed.  If type
 is set any other value then both types of loss function are evaluated.

 Since the input matrix, mat, can possess intermittent missing data,
 samonevalIM imputes these values using a logistic model.  This is the
 model specified by inmodel.

 inmodel is an NT by 6 matrix and contains a row for each time-point
 represented in the mat matrix.  Since there are no intermittent missing
 data at time 1 and at time NT these rows should be set to zero.

 A separate logistic model is fit for each timepoint t0, with t0 ranging
 from 2 to NT-1.  The dependent variable in these models is an indicator
 indicating that an individual provided data at time t0.  Only
 individuals who are on study at time t0+1 are considered in these models. 
 Setting the ith column in inmodel to 1 indicates that the following
 terms should be included in the model:

 \describe{
  \item{1}{the intercept.}
  \item{2}{the previous value of Y, that is the value of Y at t0-1.}
  \item{3}{the next available value of Y after time t0.}
  \item{4}{if the next available value of Y is at time t1, then this is
 (t1 - t0 - 1).}
  \item{5}{an interaction between the variable in column 3 and column 4.}
  \item{6}{0 if an individual is last observed at time NT, 1 if they are
  last observed at time NT - 1,  2 if they are last observed at
  time NT - 2, etc.}
 }

 For the model that runs at the penultimate time, NT-1, including all six
 terms in the model will lead to over specification.

 To include as many terms as possible in each model set inmodel to be

 \tabular{rrrrrr}{ 
 0 \tab   0  \tab  0  \tab  0  \tab  0  \tab  0\cr
 1 \tab   1  \tab  1  \tab  1  \tab  1  \tab  1\cr
 1 \tab   1  \tab  1  \tab  1  \tab  1  \tab  1\cr
 . \tab . \tab . \cr
 . \tab . \tab . \cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
}

 At a given timepoint t0, there may not be enough intermittent missing
 data to support all six terms.  In such cases variables are removed
 from the model until the regression converges.  If less than 4
 observations have intermittent missing data at a given timepoint a
 model fitting the mean is used.

 samonevalIM sequentially imputes intermittent missing data using these
 models. Once the intermittent data is imputed a loss function is
 evaluated on the data.
}
\value{
 Returns an N by 2 matrix of results if type is "H" or "F" or an N by 4
 matrix if type is "both". There is one row for each value of sigma. If
 type is "P" then the columns are sigmaH and lossH, if type is "F" then
 the columns are sigmaF and lossF, and, if type = "both" then the
 columns are sigmaH, lossH, sigmaF, and lossF.
}
\examples{
  data("VAS1")

  ## inputation moddel.
  NT                <- ncol(VAS1)
  inmodel           <- matrix(1,NT,6)
  inmodel[1,]       <- 0
  inmodel[NT,]      <- 0
  inmodel[NT-1,4:6] <- 0

  ## outcome loss function
  F1 <- samonevalIM( mat = VAS1, Npart = 10,
           sigmaList = seq(1,10,by=1),
           inmodel = inmodel,
           seed = 26,
           type = "F" )
}
