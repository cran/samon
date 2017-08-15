// ----------------------------------------------------------------------------
// Some matrix functions to support logistic regression
// Inverts positive definite matrix. 
// ----------------------------------------------------------------------------
#ifndef _samon_logit_h_
#define _samon_logit_h_

typedef struct LogisticS {
  int NT;           // number of time points (actual)
  int N0;           // number of obs in Y

  int nb;           // Max number of betas including intercept
  int Nconstr;      // number of constraints placed on betas 

  double **Data;    // data on which to perform sequential
                    // logistic regression.  Data is N0 by NT.

  int    **Last;    // last time for each subject in Data (N0 by 1)

  int    **LIter;   // Number of iterations performed in reg 
  int    **Lconvg;  // Convergence status from reg
  int    **Lnbeta;  // Number of estimates used at each t
  double **LEsts;   // Estimates
  int    **Tconstr; // constraints (0/1 if 1 beta is zero )
  int    **Tmodel;  // as Tconstr but defines the model nb by 1 
  double **LProb;   // probability of missing at each time -- used if
                    // logistic regression does not converge.

  int    **Models;  // converged models: NT by nb

  // convergence criterion
  int    MaxIter;
  double betaeps;
  double betaReleps;

  // scratch area for logistic regressions
  double **X;       // X values (N0 by nb )
  double **Y;       // Y values (0/1) (N0 by 1)
  double **Preds;   // predicted values
  int **pos;        // position in Data of X[i] and Y[i]

  double **beta0;   // input betas (nb by 1)
  double **beta1;   // updated betas (nb by 1)
  double **betau;   // update for betas (nb by 1)

  double **LT;      // lower triangle nb by nb
  double **LTi;     // lower triangle inverse
  double **D;       // first derivative of log likelyhood by beta (nb by 1)
  double **H;       // negative second derivatives (nb by nb)
  double **Hi;      // negative second derivatives inverted (nb by nb)

  double **want;    // intermittent values wanted N0 by 2
  double **have;    // predicted values N0 by 12
} LogisticS;

// make and destroy a LogisticS object  
LogisticS *initLogisticS( double **inData, int N0, int NT, int nb, int MaxIter, double betaeps, double betaReleps );
int distructLogisticS( LogisticS *logS );

// perform a logistic regression 
int samonLogit( int n, int nv, int *convg, int *Iter, double **beta0, double **beta1, LogisticS *logS );
// predicted values
int logisticP( double **X, int nr, int nc, double **beta, double **Pv);
// derivatives of the liklihood
int lkd2( double **Y, double **X, int nr, int nc, double **beta, double **lkdd, double **lkd, int **Tmodel );
#endif
