#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "basic_MatUtil.h"
#include "samon_logit.h"
#include "samon_MatUtil2.h"
#include "samon_FillRemove.h"
#include "Gen_fun.h"

void Ctime ( int t0, int *count0, int *count1, LogisticS *logS  )
{
  int i, j, cnt;
  int N0, NT, **L, **pos;
  double **Mat, **y, **x;

  N0    = logS->N0; 
  NT    = logS->NT; 
  y     = logS->Y;
  x     = logS->X;
  L     = logS->Last;
  pos   = logS->pos;
  Mat   = logS->Data;

  // --------------------------------------------------------------

  *count0 = 0;
  *count1 = 0;
  cnt = 0;
  for ( i = 0; i < N0; i++ ) { 
    if ( t0 < L[i][0] ) {
      if ( isnan(Mat[i][t0]) ) {
	y[cnt][0] = 1;
	*count1 = *count1 + 1;
      }
      else {
	y[cnt][0] = 0;
	*count0 = *count0 + 1;
      }
      x[cnt][0] = 1;
      x[cnt][1] = Mat[i][t0-1];
      x[cnt][2] = 0;
      x[cnt][3] = 0;
      x[cnt][4] = 0;
      if (( NT - 1 - L[i][0] ) > 0 ) x[cnt][5] = NT - 1 - L[i][0];
      else x[cnt][5] = 0;
      for ( j = L[i][0]; j > t0; j-- ) {
	if (!isnan(Mat[i][j])) {
	  x[cnt][2] = Mat[i][j]; 
	  if (( j - t0 ) > 1 ) {
	    x[cnt][3] = j - t0 - 1;
	    x[cnt][4] = (j - t0 - 1) * x[cnt][2];
	  }
	  else {
	    x[cnt][3] = 0;
	    x[cnt][4] = 0;
	  }
	}
      }
      pos[cnt][0] = i;
      cnt++;
    }
  }
  return;
}

int mkIM( LogisticS *logS )
{
  int i, j, N0, NT, t0;
  int **Last;
  int **Lconvg;
  double **LEsts, **Data, **LProb;
  double lp, pv, xtest, nextX, inter, nnext;

  N0      = logS->N0;
  NT      = logS->NT;
  Data    = logS->Data;
  Last    = logS->Last;
  LEsts   = logS->LEsts;
  Lconvg  = logS->Lconvg;
  LProb   = logS->LProb;

  nextX = 0.0;
  inter = 0.0;
  nnext = 0.0; 
  for ( t0 = NT - 2; t0 > 0; t0-- ) {
    if ( Lconvg[t0][0] == 1 ) {
      // the logistic converged
      for ( i = 0; i < N0; i++ ) {
	if ( t0 < Last[i][0] ) {
	  for ( j = Last[i][0]; j > t0; j-- ) { 
	    if ( !isnan(Data[i][j]) ) {
	      nextX = Data[i][j];
	      if ( j > t0 + 1 ) {
		nnext = j - t0 - 1;
		inter = (j - t0 - 1) *  Data[i][j];
	      }
	      else {
		nnext = 0;
		inter = 0;
	      }
	    }
	  }
	  lp = LEsts[t0][0] + LEsts[t0][1] * Data[i][t0-1] + LEsts[t0][2] * nextX + LEsts[t0][3] * nnext + LEsts[t0][4] * inter;
	  if (( NT - 1 - Last[i][0] ) > 0 ) lp = lp + LEsts[t0][5] * ( NT - 1 - Last[i][0]);
	  pv = 1.0 / ( 1.0 + exp(-lp) );
	  xtest = sgen();
	  if ( xtest < pv ) Data[i][t0] = NAN;
	}
      }
    }
    else {
      // else use the mean prob of intermittent missing
      for ( i = 0; i < N0; i++ ) {
	if ( t0 < Last[i][0] ) {
	  pv = LProb[t0][0];
	  xtest = sgen();
	  if ( xtest < pv ) Data[i][t0] = NAN;
	}
      }
    }
  }
  return 0;
}

int FillIn( int nv, LogisticS *logS )
{
  int t0, i, j, k, count0, count1, Total, rc;
  int NT;
  int convg, Iter, MaxIter;
  int m, c1,c2, ipick;
  int newmodel, vcnt, bcnt;
  int modelsize;
  int mcount;

  double p0, sq;

  int **DatPos;
  double **Y;
  int **Tmodel, **Tconstr;
  int Nconstr;
  double **Data, **Preds, **LEsts;
  double **beta0, **beta1, **LProb;
  double **want, **have;
  int **LIter, **Lconvg, **Lnbeta, **Models;

  NT      = logS->NT;
  Data    = logS->Data;

  Y       = logS->Y;
  Preds   = logS->Preds;

  LEsts   = logS->LEsts;
  Tmodel  = logS->Tmodel;
  Tconstr = logS->Tconstr; 
  Nconstr = logS->Nconstr;
  Models  = logS->Models;

  LIter   = logS->LIter;
  Lconvg  = logS->Lconvg;
  LProb   = logS->LProb;
  Lnbeta  = logS->Lnbeta;
 
  DatPos  = logS->pos;

  beta0   = logS->beta0;
  beta1   = logS->beta1;

  want    = logS->want;
  have    = logS->have;

  MaxIter = logS->MaxIter;

  // First let's zero these things:
  for ( t0 = 0; t0 < NT; t0++ ) {
    LProb [t0][0] = 0.0;
    Lconvg[t0][0] = 0;
    LIter [t0][0] = 0;
    Lnbeta[t0][0] = 0;
    for ( j = 0; j < nv; j++ ) {
      LEsts[t0][j] = 0.0;
    }
  }

  for ( t0 = 1; t0 < NT-1; t0++ ) {
    Ctime ( t0, &count0, &count1, logS  );
    Total = count0 + count1;

    p0 = (double) count1 / Total;

    LProb[t0][0]  = p0;
    Lconvg[t0][0] =  0;
    convg         =  0;
    LIter[t0][0]  =  0;
    Iter          =  0;

    for ( i = 0; i < nv; i++ ) beta0[i][0] = 0.0;

    if ( count1 > 3 && count0 > 3 && p0 > 0.0 && p0 < 1.0 ) {

      // base model is read from logS->Models
      Nconstr     = 0;
      modelsize   = 0;
      for ( i = 0; i < nv; i++ ) {
	if ( (logS->Models)[t0][i] == 1 ) {
	  Tmodel[i][0]  = 1;
	  Tconstr[i][0] = 0;
	  modelsize++;
	}
	else {
	  Tmodel[i][0]  = 0;
	  Tconstr[i][0] = 1;
	  Nconstr++;
	}
      }

      if ( modelsize > 0 ) {

	for ( i = 0; i < nv; i++ ) beta0[i][0] = 0.0;
	beta0[0][0] = log(p0) - log(1-p0);

	rc = samonLogit( Total, 6, &convg, &Iter, beta0, beta1, logS );

	if ( rc > 0 || convg == 0 || Iter == MaxIter ) {

	  newmodel = 0;
	  if ( rc > 0 || convg == 0 ) {
	    // try to add more constraints ( but not to the intercept )
	    for ( i = 1; i < nv; i++ ) {
	      if ( Tconstr[i][0] > 0 || fabs(beta1[i][0]) > 25.0  ) {
		if ( Tmodel[i][0] == 1 ) newmodel++; 
		Tmodel[i][0]  = 0;
		Tconstr[i][0] = 0;
	      }
	    }
	  }

	  if ( (rc > 0 || convg == 0 || Iter == MaxIter ) && newmodel == 0 ) {
	    // didn't work -- could be MaxIter or the intercept is wrong 
	    // dump the last variable
	    vcnt = 0;
	    for (i = 0; i < nv; i++ ) {
	      if ( Tmodel[i][0] == 1 ) vcnt = i;
	    }
	    if ( vcnt > 0 ) {
	      Tmodel[vcnt][0]  = 0;
	      Tconstr[vcnt][0] = 0;
	    }
	  }

	  do {
	    modelsize--;

	    for ( i = 0; i < nv; i++ ) beta0[i][0] = 0.0;
	    beta0[0][0]   = log(p0) - log(1-p0);
	    Lconvg[t0][0] = 0;
	    convg         = 0;
	    LIter[t0][0]  = 0;
	    Iter          = 0;

	    rc = samonLogit( Total, 6, &convg, &Iter, beta0, beta1, logS );

	    newmodel = 0;
	    if ( rc > 0 || convg == 0 ) {
	      // try and reduce the model: by adding more constraints (again not the int)
	      for ( i = 1; i < nv; i++ ) {
		if ( Tconstr[i][0] > 0 || fabs(beta1[i][0]) > 25.0  ) {
		  if ( Tmodel[i][0] == 1 ) newmodel++; 
		  Tmodel[i][0]  = 0;
		  Tconstr[i][0] = 0;
		}
	      }
	    }
	    if ( (rc > 0 || convg == 0) && newmodel == 0 ) {
	      // didn't work -- dump the last variable
	      vcnt = 0;
	      for (i = 0; i < nv; i++ ) {
		if ( Tmodel[i][0] == 1 ) vcnt = i;
	      }
	      if ( vcnt > 0 ) {
		Tmodel[vcnt][0]  = 0;
		Tconstr[vcnt][0] = 0;
	      }
	    }

	  } while ( (rc > 0 || convg == 0 || Iter == MaxIter) && modelsize > 0 ); 
	  modelsize++;
	}
      }

      bcnt = 0;
      for ( i = 0; i < nv; i++ ) {
	if ( Tmodel[i][0] == 1 ) {
	  bcnt++;
	  Models[t0][i] = 1;
	}
	else {
	  Models[t0][i] = 0;
	}
      } 

      Lconvg[t0][0] = convg;
      LIter[t0][0]  = Iter;
      Lnbeta[t0][0] = bcnt;
      if ( convg == 1 ) {
	for ( j = 0; j < nv; j++ )  {
	  LEsts[t0][j]  = beta1[j][0];
	}
      }
      else {
	for ( j = 0; j < nv; j++ ) 
	  LEsts[t0][j] = 0.0;
      }
    }

    if ( convg == 0 ) {
      // No logistic done or it didn't converge
      bcnt = 0;
      for ( i = 0; i < nv; i++ ) {
	Models[t0][i] = 0;
      } 
    }

    if ( count1 == 0 ) {
      // nothing to do
    }
    else if ( convg == 0 ) {
      if ( p0 == 1.0 ) {
	// problem -- do nothing?
      }
      else {
	// randomly pick values 
	c1 = 0;
	c2 = 0;
	for ( i = 0; i < Total; i++ ){

	  if ( Y[i][0] == 1 ) {
	    want[c1][1] = DatPos[i][0];
	    c1++;
	  }
	  else {
	    have[c2][1] = DatPos[i][0];
	    c2++;
	  }
	}

	for ( i = 0; i < c1; i ++ ) {
	  ipick = (int) floor(c2 * sgen());
	  if ( ipick > (c2 - 1)) ipick = c2 - 1;
	  // Note change here i ==> ipick
	  Data[ (int) want[i][1] ][t0] = Data[ (int) have[ipick][1] ][t0];
	}
      }
    }
    else {

      c1 = 0;
      c2 = 0;
      for ( i = 0; i < Total; i++ ){

	if ( Y[i][0] == 1 ) {
	  want[c1][0] = Preds[i][0];
	  want[c1][1] = DatPos[i][0];
	  for ( j = 2; j < 12; j=j+2 ) {
	    want[c1][j  ] = 1E50;
	    //	    want[c1][j+1] = R_NaReal;
	    want[c1][j+1] = NAN;
	  }
	  c1++;
	}
	else {
	  have[c2][0] = Preds[i][0];
	  have[c2][1] = DatPos[i][0];
	  c2++;
	}
      }

      for ( i = 0; i < c1; i++ ) {
	for ( j = 0; j < c2; j++ ) {
	  sq = pow(want[i][0] - have[j][0],2.0);
	  for ( m = 2; m < 12; m = m + 2 ) {
	    if ( sq < want[i][m] ) {
	      for ( k = 10; k > m; k=k-2 ) {
		want[i][k]   = want[i][k-2];
		want[i][k+1] = want[i][k-1];
	      }
	      want[i][m]   = sq;
	      want[i][m+1] = have[j][1];
	      m = 12;
	    }
	  }
	}
      }

      for ( i = 0; i < c1; i++ ) {

	// how many to pick from:
	mcount = 0;
	for ( m = 2; m < 12; m = m + 2 ) {
	  if ( !isnan(want[i][m]) && !isnan(want[i][m+1]) ) mcount++;
	}

	ipick = (int) floor(mcount * sgen());
	if ( ipick > mcount ) ipick = mcount;
	ipick = 3 + 2 * ipick;

	Data[ (int) want[i][1] ][t0] = Data[ (int) want[i][ipick] ][t0];
      }

    }
  }
  return nv;
}
