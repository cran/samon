\name{samonGen}
\alias{samonGen}
\title{Sensitivity Analysis for monotone missing data}
\description{Given data from one arm of a repeated measures clinical
  trial, produces a sample based on the optimal smoothing parameters.}
\usage{samonGen(mat, Npart = 10, InitialSigmaH = 1.0, HighSigmaH = 2.0,
InitialSigmaF = 1.0, HighSigmaF = 2.0, NSamples = 0, seed0 = 1,
MaxIter = 25, FAconvg = 1E-7, FRconvg = 1E-7, SAconvg = 1E-7)
}
\arguments{
  \item{mat}{ a matrix with the (i,j) entry representing the outcome
    value for subject i at time-point j.}
  \item{Npart}{ Number of partitions to use when estimating optimal
    smoothing parameters, sigma H and sigma F.}
  \item{InitialSigmaH}{ Initial value when calculating optimal sigma H. }
  \item{HighSigmaH}{ Upper bound of search region when calculating optimal sigma H. }
  \item{InitialSigmaF}{ Initial value when calculating optimal sigma F. }
  \item{HighSigmaF}{ Upper bound of search region when calculating optimal sigma F. }
  \item{NSamples}{ Number of parametric bootstrap samples to generate. }
  \item{seed0}{ Seed to use. }
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
 FAconvg:  stop when abs( fsub(i+1) - f_i ) < FAconvg, where f_i is
 the loss function value at iteration i.
 FRconvg: stop when abs( (fsub(i+1) - f_i)/(fsub(i+1) + f_i)) < FRconvg.
 SAconvg: stop when the absolute step size falls below SAconvg,
   i.e., abs( xsub(i+1) - x_i ) < SAconvg.
 HighSigmaH:  should the value of sigma H go above this value,
 then the optimal value of sigma H is set to HighSigmaH.  This is useful
 if larger values of sigma H do not change the missingness model
 substantially.
 HighSigmaF: should the value of sigma F go above this value, then
 the optimal value of sigma F is set to HighSigmaF.  This is useful
 if larger values of sigma F do not change the missingness model
 substantially.
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
data(samonPANSS1)

Sample1 <- samonGen(
  mat            = samonPANSS1,
  Npart          =           5, # number of partitions
    
  InitialSigmaH  =        15.6, # initial value
  HighSigmaH     =        25.0, # high value for H
    
  InitialSigmaF  =         8.6, # initial value
  HighSigmaF     =        15.0, # high value for F
    
  seed           =         211,
  NSamples       =          30 )
}
