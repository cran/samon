\name{samon}
\alias{samon}
\title{Sensitivity Analysis for monotone missing data}
\description{Given data from one arm of a repeated measures clinical
  trial, produces estimates of the expected value of the outcome at the
  final time-point for a range of sensitivity parameters.}
\usage{samon(mat, Npart = 10, InitialSigmaH = 1.0, HighSigmaH = 2.0,
InitialSigmaF = 1.0, HighSigmaF = 2.0, lb = 0, ub = 101, 
zeta1 = 1, zeta2 = 1,NSamples = 0, seed0 = 1, MaxIter = 25,
FAconvg = 1E-7, FRconvg = 1E-7, SAconvg = 1E-7, alphaList= c(0), 
MJackknife = FALSE, SJackknife = FALSE, retIFiM = FALSE, 
retIFiS = FALSE, Tfun= NULL)
}
\arguments{
  \item{mat}{ a matrix with the (i,j) entry representing the outcome
    value for subject i at time-point j.} 
  \item{Npart}{ Number of partitions to use when estimating optimal
    smoothing parameters, sigma H and sigma F.}
  \item{InitialSigmaH}{ Initial value when calculating optimal sigma H. }
  \item{HighSigmaH}{ Upper bound of search region when calculating optimal sigma H. }
  \item{InitialSigmaF}{ Initial value when calculation optimal sigma F. }
  \item{HighSigmaF}{ Upper bound of search region when calculating optimal sigma F. }
  \item{lb}{ Lower bound for Y. }
  \item{ub}{ Upper bound for Y. }
  \item{zeta1}{ parameter to cumulative beta. }
  \item{zeta2}{ parameter to cumulative beta. }
  \item{NSamples}{ Number of parametric bootstrap samples to generate. }
  \item{seed0}{ Seed to use. }
  \item{MaxIter}{ Maximum iterations to use in optimizer. }
  \item{FAconvg}{ Absolute change in function convergence criterion. }
  \item{FRconvg}{ Relative change in function convergence criterion. }
  \item{SAconvg}{ Step size convergence criterion. }
  \item{alphaList}{ a vector of sensitivity parameters }
  \item{MJackknife}{ Jackknife main data (logical)}
  \item{SJackknife}{ Jackknife bootstrap samples (logical) }
  \item{retIFiM}{ return individual IF estimates from main data. }
  \item{retIFiS}{ return individual IF estimates from samples. }
  \item{Tfun}{ n by 2 matrix with tilting function values. }
}

\details{
 The matrix mat represents repeated measure outcome data from a single
 arm or treatment group of a trial.  Each row represents the data from a
 single subject and each column data from a single time-point.

 The values in the first column of mat are the baseline values and should
 not have missing values.  Samon creates one-step bias corrected estimates of
 the mean value of Y at the last time-point for a number of specified
 sensitivity parameters, alpha.

 Samon determines two smoothing parameters, sigma H, which represents
 smoothing in the "missingness" model and, sigma F, which represents
 smoothing in the "outcome" model.  These smoothing parameters are
 determined by minimizing loss functions.  Minimization is performed
 using Newton's method. The parameter InitialSigmaH is used as the
 initial value in the optimization of the missingness model.  In a
 similar vein InitialSigmaF is used as the initial value in the
 optimization of the outcome model.

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
 
 Using the optimal values for sigma H and sigma F, samon produces
 estimates of the mean value of the outcome at the last visit, for the
 set of sensitivity parameters given in the vector alphaList.

 If the sensitivity parameters are applied to the cumulative beta function,
 IB.  Four parameters are used to map the values of the outcome, y to the
 domain of the  cumulative beta. r(y) = IB( (y-lb)/(ub-lb), zeta1,
 zeta2) where zeta1 and zeta2 are the shape parameters.

 Alternatively samon may be called with a user specified sensitivity
 function. This is given as an n by 2 matrix with the first column
 holding the unique values of the data in mat and the second column the
 value of the sensitivity function at that value.
}
\value{
 samon returns a list which includes the following:
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

  \item{IFM}{ matrix with 7 columns containing one-step bias-corrected
    estimates (sometimes referred to here as IF estimates). Columns are
    sample, type, alpha, est0, se0, IFEst, IFVar.  The columns sample
    and type take the value 0.  IFEst and IFVar are the one-step bias
    corrected estimates and an estimate of their variance for the
    original data, mat, at each value of the sensitivity parameters, alpha. }

  In a similar manner, "HS", "FS", and "IFS" contain results from
  parametric bootstraps (should NSamples be greater than 0).  The first
  columns in this case contain the sample number, 1, 2, ...,
  NSamples. In other respects these matrices have the same structure as
  the corresponding matrix "HM", "FM" and "IFM".

  If MJackknife is set to TRUE, then, "HMjk", "FMjk", "IFMjk" are
  produced, holding the jackknife results for the main data mat.  Again
  the first column "sample" is set to 0, and the second column is termed
  "Dropped" and gives the observation number dropped to form the
  jackknife sample.

  Finally, if SJackknife is set to TRUE, then, "HSjk", "FSjk", "IFSjk"
  are produced.  These hold the jackknife results for the bootstrap
  samples.
}
\seealso{
  The samon_userDoc.pdf file in the Examples subdirectory.
}  
\examples{
data("samonPANSS1")

Results1 <- samon(
  mat            = samonPANSS1,
  Npart          =          10, # number of partitions
    
  InitialSigmaH  =        15.0, # initial value
  HighSigmaH     =        50.0, # high value for H
    
  InitialSigmaF  =         8.0, # initial value
  HighSigmaF     =        50.0, # high value for F
    
  lb             =          30, # parameters to
  ub             =         210, # cumulative 
  zeta1          =         4.0, # beta distribution
  zeta2          =         7.0,

  alphaList      =        -1:1 )
}
