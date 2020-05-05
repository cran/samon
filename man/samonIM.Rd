\name{samonIM}
\alias{samonIM}
\title{samonIM: Sensitivity analysis for monotone missing and
  intermittent missing data}
\description{Given data from one arm of a repeated measures clinical
  trial, produces estimates of the expected value of the outcome at the
  final time-point for a range of sensitivity parameters.}
\usage{samonIM(mat, Npart = 10, InitialSigmaH = 1.0, 
HighSigmaH = 2.0, InitialSigmaF = 1.0, HighSigmaF = 2.0, 
inmodel = inmodel, NSamples = 0, NIMimpute = 1, lb = 0, 
ub = 101, zeta1 = 1, zeta2 = 1, seed0 = 1, seed1 = 1, 
MaxIter = 25, FAconvg = 1E-7, FRconvg = 1E-7, 
SAconvg = 1E-7, alphaList = c(0), retIFiM  = FALSE, 
retIFiS  = FALSE, retSample = FALSE,
retFMatM = FALSE, retFMatS = FALSE, Tfun=NULL)}
\arguments{
  \item{mat}{ a matrix with the (i,j) entry representing the outcome
    value for subject i at time-point j. }
  \item{Npart}{ Number of partitions to use when estimating optimal
    smoothing parameters, sigma H and sigma F.}
  \item{InitialSigmaH}{ Initial value when calculating optimal sigma H. }
  \item{HighSigmaH}{ Upper bound of search region when calculating optimal sigma H. }
  \item{InitialSigmaF}{ Initial value when calculation optimal sigma F. }
  \item{HighSigmaF}{ Upper bound of search region when calculation
    optimal sigma F. }
  \item{NSamples}{ Number of parametric bootstrap samples to generate. }
  \item{inmodel}{ NT by 6 matrix indicating which variables to include
    in logistic model. }
  \item{NIMimpute}{ How many times to fill in intermittent missing data. } 
  \item{lb}{ Lower bound for Y. }
  \item{ub}{ Upper bound for Y. }
  \item{zeta1}{ parameter to cumulative beta. }
  \item{zeta2}{ parameter to cumulative beta. }
  \item{seed0}{ Seed to use for parametric bootstrap sampling }
  \item{seed1}{ Seed to use for filling intermittent missing data }
  \item{MaxIter}{ Maximum iterations to use in optimizer. }
  \item{FAconvg}{ Absolute change in function convergence criterion. }
  \item{FRconvg}{ Relative change in function convergence criterion. }
  \item{SAconvg}{ Step size convergence criterion. }
  \item{alphaList}{ a vector of sensitivity parameters }
  \item{retIFiM}{ return individual IF estimates from main data. }
  \item{retIFiS}{ return individual IF estimates from sample data. }
  \item{retSample}{ return the Sample generated }
  \item{retFMatM}{ return the main data with intermittent missing data
    filled in NIMimpute + 1 times}
  \item{retFMatS}{ retrun the sample data with intermittent missing data
    filled in NIMimpute + 1 times}
  \item{Tfun}{ n by 2 matrix with tilting function values. }
}

\details{
 The matrix mat represents repeated measure outcome data from a single
 arm or treatment group of a trial.  Each row represents the data from a
 single subject and each column data from a single time-point.  

 The values in the first column of mat are the baseline values and should
 not be missing.  Samon creates bias-corrected estimates of the mean
 value of Y at the last time-point for a number of specified sensitivity
 parameters, alpha.

 Samon determines two smoothing parameters, sigma H, which represents
 smoothing in the "missingness" model and, sigma F, which represents
 smoothing in the "outcome" model.  These smoothing parameters are
 determined by minimizing loss functions.  Minimization is performed
 using Newton's method. The parameter InitialSigmaH is used as the
 initial value in the optimization of the missingness model and
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

 The inmodel matrix specifies a set of logistic models used to impute
 intermittent missing data.  This is an NT by 6 matrix of 0s and 1s.
 Each row represents a time-point.  Since intermittent missing data
 cannot occur at baseline or at the final timepoint, the first and NT'th
 row of inmodel are 0.  The six values on a row indicate the following
 terms should be included in the model:
 \describe{
 \item{1}{ a value of 1 in this position indicates that the logistic
 model should include an intercept}
 \item{2}{ a value of 1 in this position indicates that the model should
 contain a term for the previous available outcome}
 \item{3}{ a value of 1 in this position indicates that the model should
 contain a term for the next available outcome}
 \item{4}{ a value of 1 in this position indicates that the model should
 contain a term for the location of the next available outcome}
 \item{5}{ a value of 1 in this position indicates that the model should
 contain a term for the interaction between the next available outcome
 and the position of the next available outcome}
 \item{6}{ a value of 1 in this position indicates that the model should
   contain a term for the last available outcome before dropout.}
 }

 For example suppose the data consists of the 8 columns v1-v8 and
 suppose we wish to model the intermittent missingness at the third
 time-point, that is, we wish to model when v3 has an intermittent
 missing value.   We might model this as is.na(v3) ~ 1 + v2 + n3 + i3 +
 i3 * n3 + l3 where n3 refers to the next non-missing value of v after
 the third time-point. i3 refers to the location of that non-missing
 value, i3 * n3 and interaction between these and l3 the last available
 value of v for an individual.  Specifically we might calculate n3, i3
 and l3 as shown:

 \tabular{rrrrrrrrrrrr}{
 v1 \tab  v2 \tab    v3 \tab    v4  \tab   v5  \tab   v6 \tab    v7  \tab   v8  \tab     v2 \tab    n3 \tab    i3 \tab    l3\cr
 10 \tab     8 \tab    NA \tab     9 \tab    32 \tab    NA \tab    NA \tab    NA \tab       8 \tab     9 \tab     1 \tab    32\cr
 11 \tab     4 \tab     3 \tab     5 \tab     4 \tab     7 \tab     7 \tab    NA \tab       4 \tab     5 \tab     1 \tab     7\cr
 22 \tab    15 \tab    NA \tab    NA \tab    NA \tab     6 \tab     2 \tab     9 \tab      15 \tab     6 \tab     3 \tab     9\cr
}
 
 Note that when considering missing values at the NT - 1 time-point the
 last three terms in the model are not identifiable.

 Fitting as many terms in the intermittent missing model might use an
 inmodel parameter with a matrix like this:

 \tabular{rrrrrr}{ 
 0 \tab   0  \tab  0  \tab  0  \tab  0  \tab  0\cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   1 \tab   1 \tab   1\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
}

 A simpler model using only the intercept and the previous and next
 available non-missing value would be specified by

 \tabular{rrrrrr}{ 
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
}

 This simpler model but restricted to an intercept only model at the
 second and fourth timepoints can be specified by using:  
 
 \tabular{rrrrrr}{ 
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
 1 \tab   0 \tab   0  \tab  0  \tab  0  \tab  0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 1 \tab   1 \tab   1 \tab   0 \tab   0 \tab   0\cr
 0 \tab   0 \tab   0 \tab   0 \tab   0 \tab   0\cr
}

 In practice the amount of intermittent missing data at a timepoint may
 not be large enough to support a logistic model with all six
 parameters.  At each step samonIM  will attempt to fit the models
 indicated by the inmodel matrix and reduce the number of model
 parameters if collinearity is found among the independent variables.  
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
\enumerate{
 \item absolute function convergence was met 
 \item relative function convergence was met 
 \item second derivative has become too small
 \item maximum iterations reached 
 \item value reset to HighSigmaH or HighSigmaF 
 \item loss function smaller at HighSigmaH or HighSigmaF 
}

  \item{Sample}{ The generated sample of dimension NSamples by ncol(Mat)}
}
\seealso{
  The samon_userDoc.pdf file in the Examples subdirectory.
}  
\examples{
data("VAS1")

# inputation model
NT            <- ncol(VAS1)
inmodel       <- matrix(1,NT,6)
inmodel[1,]   <- 0
inmodel[NT,]  <- 0
inmodel[,4:6] <- 0

Results <- samonIM(
  mat           =    VAS1, # imput matrix
  Npart         =       2, # number of partitions
    
  InitialSigmaH =    25.0, # initial value
  HighSigmaH    =   100.0, # high value for H
    
  InitialSigmaF =     8.0, # initial value
  HighSigmaF    =   100.0, # high value for F
    
  lb            =       0, # parameters for
  ub            =     102, # cumulative 
  zeta1         =     1.2, # beta distribution
  zeta2         =     1.6,

  NSamples      =       0, # no of bootstraps
  NIMimpute     =       2, # no of imputations 
    
  seed0         =     441, # seed for bootstraps
  seed1         =     511, # seed for imputations
  inmodel       = inmodel, # input model
    
  alphaList     =   -1:1 )
}
