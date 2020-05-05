\name{samonGenIM}
\alias{samonGenIM}
\title{Sensitivity Analysis for monotone missing data}
\description{Given data from one arm of a repeated measures clinical
  trial with intermittent missing data, produces a sample based on the
  optimal smoothing parameters.}
\usage{samonGenIM(mat, Npart = 10, InitialSigmaH = 1.0, 
HighSigmaH = 2.0, InitialSigmaF = 1.0, HighSigmaF = 2.0, 
inmodel = inmodel, NSamples = 0, seed0 = 1, seed1 = 1, 
MaxIter = 25, FAconvg = 1E-7, FRconvg = 1E-7,SAconvg = 1E-7)
}
\arguments{
  \item{mat}{ a matrix with the (i,j) entry representing the outcome
    value for subject i at time-point j.}
  \item{Npart}{ Number of partitions to use when estimating optimal
    smoothing parameters, sigma H and sigma F.}
  \item{InitialSigmaH}{ Initial value when calculating optimal sigma H. }
  \item{HighSigmaH}{ Upper bound of search region when calculating optimal sigma H. }
  \item{InitialSigmaF}{ Initial value when calculating optimal sigma F. }
  \item{HighSigmaF}{ Upper bound of search region when calculating
    optimal sigmaF. }
  \item{inmodel}{ NT by 6 matrix indicating which variables to include
    in logistic model. }
  \item{NSamples}{ Number of parametric bootstrap samples to generate. }
  \item{seed0}{ Seed to use when generating bootstrap samples. }
  \item{seed1}{ Seed to use when imputing intermittent missing data. }
  \item{MaxIter}{ Maximum iterations to use in optimizer. }
  \item{FAconvg}{ Absolute change in function convergence criterion. }
  \item{FRconvg}{ Relative change in function convergence criterion. }
  \item{SAconvg}{ Step size convergence criterion. }
}

\details{
 The matrix mat represents repeated measure outcome data from a single
 arm or treatment group of a trial.  Each row represents the data from a
 single subject and each column data from a single time-point. 

 The values in the first column of mat are the baseline values and should
 not be missing.
 
 Samon determines two smoothing parameters, sigma H, which represents
 smoothing in the "missingness" model and, sigma F, which represents
 smoothing in the "outcome" model.  These smoothing parameters are
 determined by minimizing loss functions.  Minimization is performed
 using Newton's method. The parameter InitialSigmaH is used as the
 initial value in the optimization of the missingness model, and,
 InitialSigmaF is used as the initial value in the optimization of the
 outcome model.

 A number of stopping criteria are available:
 MaxIter:  the maximum number of iterations to perform 
 FAconvg:  stop when abs( f_{i+1} - f_i ) < FAconvg, where f_i is
 the loss function value at iteration i.
 FRconvg: stop when abs( (f_{i+1} - f_i)/(f_{i+1} + f_i) ) < FRconvg.
 SAconvg: stop when the absolute step size falls below SAconvg,
   i.e., abs( x_{i+1} - x_i ) < SAconvg.
 HighSigmaH:  should the value of sigma H go above this value,
 then the optimal value of sigma H is set to HighSigmaH.  This is useful
 if larger values of sigma H do not change the missingness model
 substantially.
 HighSigmaF: should the value of sigma F go above this value, then
 the optimal value of sigma F is set to HighSigmaF.  This is useful
 if larger values of sigma F do not change the missingness model
 substantially.


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

 samonGenIM sequentially imputes intermittent missing data using these
 models. Once the intermittent data is imputed a loss function is
 evaluated on the data.
}
\value{
 samonGen returns a list which includes the following:
  \item{HM}{ matrix of results from sigma H optimization for the input
 data mat. Columns are sample, type, return code, iterations, optimal
 sigma H, and the loss function value at optimal sigma H. The "M" in the
 name "HM" refers to the main or input matrix mat.  In this case the
 sample and type columns are set to 0.}
 
  \item{FM}{ matrix of results from sigma F optimization for the input
 data mat. Columns are sample, type, return code, iterations, optimal
 sigma F, and loss function value at optimal sigma F. The "M" in the
 name "FM" refers to the main or input matrix mat.  In this case the
 sample and type columns are set to 0.}

 The return code takes the following values:   
  \item{1}{ absolute function convergence was met }
  \item{2}{ relative function convergence was met }
  \item{3}{ second derivative has become too small }
  \item{4}{ maximum iterations reached }
  \item{5}{ value reset to HighSigmaH or HighSigmaF }
  \item{6}{ loss function smaller at HighSigmaH or HighSigmaF }

 \item{Sample}{ a matrix of size Nsamples by NT with the sample }
}
\seealso{
  The samon_userDoc.pdf file in the Examples subdirectory.
}  
\examples{
data("VAS1")

# define the imputation models
NT                <- ncol(VAS1)
inmodel           <- matrix(1,NT,6)
inmodel[1,]       <- 0
inmodel[NT,]      <- 0
inmodel[NT-1,4:6] <- 0

Sample1 <- samonGenIM(
  mat           =    VAS1, # imput matrix
  Npart         =       5, # number of partitions
    
  InitialSigmaH =    27.9, # initial value
  HighSigmaH    =   100.0, # high value for H
    
  InitialSigmaF =     7.3, # initial value
  HighSigmaF    =   100.0, # high value for F
    
  NSamples      =      12, # number of samples
    
  seed0         =     441, # sample seed
  seed1         =     511, # impute seed 
  inmodel       = inmodel) # input model
}
