#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "samon_logit.h"
#include "basic_MatUtil.h"
#include "samon_MatUtil2.h"

// make a LogisticS object
LogisticS *initLogisticS( double **inData, int N0, int NT, int nb, int MaxIter, double betaeps, double betaReleps )
{
  int i,j;
  LogisticS *logS;

  logS = malloc( sizeof( LogisticS) );

  logS->N0         = N0;
  logS->NT         = NT;
  logS->nb         = nb;

  logS->Data       = mkMatd(N0,NT);
  logS->Last       = mkMati(N0, 1);

  logS->LEsts      = mkMatd(NT,nb);
  logS->Models     = mkMati(NT,nb);
  logS->Tmodel     = mkMati(nb, 1);
  logS->Tconstr    = mkMati(nb, 1);
  logS->LProb      = mkMatd(NT, 1);
  logS->LIter      = mkMati(NT, 1);
  logS->Lconvg     = mkMati(NT, 1);
  logS->Lnbeta     = mkMati(NT, 1);

  logS->MaxIter    = MaxIter;
  logS->betaeps    = betaeps;
  logS->betaReleps = betaReleps;

  logS->X          = mkMatd(N0,nb);
  logS->Y          = mkMatd(N0, 1);
  logS->Preds      = mkMatd(N0, 1);
  logS->pos        = mkMati(N0, 1);

  logS->beta0      = mkMatd(nb, 1);
  logS->beta1      = mkMatd(nb, 1);
  logS->betau      = mkMatd(nb, 1);

  logS->LT         = mkMatd(nb,nb);
  logS->LTi        = mkMatd(nb,nb);
  logS->D          = mkMatd(nb, 1);
  logS->H          = mkMatd(nb,nb);
  logS->Hi         = mkMatd(nb,nb);

  logS->want       = mkMatd(N0,12);
  logS->have       = mkMatd(N0, 2);

  for ( i = 0; i < N0; i++ ) {
    for ( j = 0; j < NT; j++ ) {
      (logS->Data)[i][j] = inData[i][j];
      if ( !isnan(inData[i][j]) ) (logS->Last)[i][0] = j;
    }
  }

  // zero these:
  for ( i = 0; i < NT; i++ ) {
    (logS->LProb)[i][0] = 0.0;
    (logS->Lconvg)[i][0] = 0;
    (logS->LIter )[i][0] = 0;
    (logS->Lnbeta)[i][0] = 0;
    for ( j = 0; j < nb; j++ ) {
      (logS->LEsts )[i][j] = 0.0;
      (logS->Models)[i][j] = 0;
    }
  }
  return logS;
}

// destroy a LogisticS object
int distructLogisticS( LogisticS *logS )
{
  freeMatd(  logS->Data    );
  freeMati(  logS->Last    );
  freeMati(  logS->LIter   );
  freeMati(  logS->Lconvg  );
  freeMati(  logS->Lnbeta  );

  freeMatd(  logS->LEsts   );
  freeMati(  logS->Models  );
  freeMati(  logS->Tmodel  );
  freeMati(  logS->Tconstr );
  freeMatd(  logS->LProb   );

  freeMatd(  logS->X       );
  freeMatd(  logS->Y       );
  freeMatd(  logS->Preds   );
  freeMati(  logS->pos     );

  freeMatd(  logS->beta0   );
  freeMatd(  logS->beta1   );
  freeMatd(  logS->betau   );

  freeMatd(  logS->LT      );
  freeMatd(  logS->LTi     );
  freeMatd(  logS->D       );
  freeMatd(  logS->H       );
  freeMatd(  logS->Hi      );

  freeMatd(  logS->want    );
  freeMatd(  logS->have    );


  free(logS);
  return 0;
}

int logisticP( double **X, int nr, int nc, double **beta, double **Pv)
{
  int i, j;

  for ( i=0; i<nr; i++ ) {
    Pv[i][0] = 0.0;  
    for ( j=0; j<nc; j++ ) {
      Pv[i][0] = Pv[i][0] + beta[j][0] * X[i][j];
    }
    // avoid under or overflow
    if ( Pv[i][0] < -700.0 ) Pv[i][0] = -700.0;
    else if ( Pv[i][0] > 700.0 ) Pv[i][0] = 700.0;
    Pv[i][0] = 1.0 / ( 1.0 + exp(-Pv[i][0]) );
  }
  return 0;
}

int lkd2( double **Y, double **X, int nr, int nc, double **beta, double **lkdd, double **lkd, int **Tmodel)
{
  int i,j,k;
  double pv;

  for ( i=0; i<nc; i++ ) {
    for ( j=0; j<nc; j++ ) lkdd[i][j] = 0.0;
    lkd[i][0] = 0.0;
  }

  for ( i=0; i<nr; i++ ) {
    pv = 0.0;  
    for ( j=0; j<nc; j++ ) {
      if ( Tmodel[j][0] == 1 ) {
	pv = pv + beta[j][0] * X[i][j];
      }
    }
    // avoid under or overflow
    if ( pv < -700.0 ) pv = -700.0;
    else if ( pv > 700.0 ) pv = 700.0;
    pv = 1.0 / ( 1.0 + exp(-pv) );
    for ( j = 0; j < nc; j++ ) {
      if ( Tmodel[j][0] == 1 ) {
	lkd[j][0]  = lkd[j][0] + ( Y[i][0] - pv ) * X[i][j];
	for ( k = 0; k < nc; k++ ) {
	  if ( Tmodel[k][0] == 1 ) {
	    lkdd[j][k] = lkdd[j][k] + (1-pv)*pv*X[i][j]*X[i][k];
	  }
	}
      }
    }
  }

  // these are not in the model
  for ( j=0; j<nc; j++ ) {
    if ( Tmodel[j][0] == 0 ) {
      for ( k=0; k < nc; k++ ) {
	lkdd[j][k] = 0.0;
	lkdd[k][j] = 0.0;
      } 
      lkdd[j][j] = 1.0;
      lkd[j][0]  = 0.0;
    }
  }
  return 0;
}

int samonLogit( int n, int nv, int *convg, int *Iter, double **beta0, double **beta1, LogisticS *logS )
{
  int i;
  int iter, rc;
  int **Tmodel, **Tconstr;
  int Npivits;
  double betaDiff, betaRelDiff, maxBeta;
  int MaxIter = 25;
  double betaeps = 10E-10;
  double betaReleps = 10E-10;

  double **H, **Hi, **D, **betaup;
  double **Y,  **X,  **Preds;
  double **LT, **LTi;

  H           = logS->H;
  Hi          = logS->Hi;
  D           = logS->D;
  betaup      = logS->betau;

  Tmodel      = logS->Tmodel;
  Tconstr     = logS->Tconstr; 

  Y           = logS->Y;
  X           = logS->X;
  Preds       = logS->Preds;

  LT          = logS->LT;
  LTi         = logS->LTi;

  MaxIter     = logS->MaxIter;
  betaeps     = logS->betaeps;
  betaReleps  = logS->betaReleps;

  // before we start set beta1 to beta0
  for ( i = 0; i < nv; i++ ) beta1[i][0] = beta0[i][0];

  // first and negative second ders of the likelihood
  rc = lkd2( Y, X, n, nv, beta0, H, D, Tmodel);
  rc = Matinv( H, nv, Hi, LT, LTi, &Npivits, Tmodel, Tconstr );

  if ( rc > 0 ) {
    // Not of full rank or other problem
    *convg = 0;
    return 1;
  }

  rc = MatMult( Hi, nv, nv, D, nv, 1, betaup);

  maxBeta     = 0.0;
  betaDiff    = 0.0;
  betaRelDiff = 0.0;
  for ( i = 0; i < nv; i++ ) {
    if ( fabs(betaup[i][0]) > betaDiff ) betaDiff = fabs(betaup[i][0]);
    beta1[i][0] = beta0[i][0] + betaup[i][0];
    if ( fabs(beta1[i][0]) > maxBeta ) maxBeta = fabs(beta1[i][0]);
  }

  // if change in the beta is too big
  if ( (betaDiff > 25.0) || (maxBeta > 25.0)) {
    *convg = 0;
    return 2;
  }

  // -----------------------

  iter = 1; 
  do {
    rc = lkd2( Y, X, n, nv, beta1, H, D, Tmodel);
    rc = Matinv( H, nv, Hi, LT, LTi, &Npivits, Tmodel, Tconstr );

    if ( rc > 0 ) {
      // Not of full rank or other problem
      *convg = 0;
      return 3;
    }
    rc = MatMult( Hi, nv, nv, D, nv, 1, betaup);

    betaDiff    = 0.0;
    betaRelDiff = 0.0;
    maxBeta     = 0.0;
    for ( i = 0; i < nv; i++ ) {
      if ( fabs(betaup[i][0]) > betaDiff ) betaDiff = fabs(betaup[i][0]);
      if ( fabs( beta1[i][0] + 0.5 * betaup[i][0] ) > 0.0 ) {
	if ( betaRelDiff < fabs(betaup[i][0]) / fabs( beta1[i][0] + 0.5 * betaup[i][0]) )
	  betaRelDiff = fabs(betaup[i][0]) / fabs( beta1[i][0] + 0.5 * betaup[i][0]);
      }
      beta1[i][0] = beta1[i][0] + betaup[i][0];
      if ( fabs(beta1[i][0]) > maxBeta ) maxBeta = fabs(beta1[i][0]);
    }

    // if change in the beta is too big
    if ( (betaDiff > 25.0) || (maxBeta > 25.0)) {
      *convg = 0;
      return 4;
    }
    iter++; 
  } while (iter < MaxIter && betaDiff > betaeps && betaRelDiff > betaReleps );
  *Iter = iter;

  rc = logisticP( X, n, nv, beta1, Preds);

  if ( iter < MaxIter ) *convg = 1;
  else *convg = 0;

  if ( *convg == 0 ) return 5;
  else return 0;
}
