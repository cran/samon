// ----------------------------------------------------------------------------
// a collection of functions to write data to output arrays
// and a function to set up the baseline logistic model
// ----------------------------------------------------------------------------

#include "samon_Util.h"

// set the baseline logistic model
int setBase(int **model, int nt, int nb ) {
  int i, j;

  for ( i = 0; i < nt; i++ ) {
    for ( j = 0; j < nb; j++ ) {
      if ( i > 0 && i < nt - 2 ) {
	model[i][j] = 1;
      }
      else if ( i == nt - 2 && j < 3 ) {
	model[i][j] = 1;
      }
      else {
	model[i][j] = 0;
      }
    }
  }
  return 0;
}

// write the logistic estimates and models to output 
int logStoOut(LogisticS *logS,double **LEptr2, int **Modelptr2, int nt, int nb, int id, int id2 ) {
  int i, j;
  double *LEstsptr;
  int    * Modelsptr;

  LEstsptr  = *LEptr2; 
  for ( i = 0; i < nt; i++ ) { 
    *LEstsptr++ = id;
    *LEstsptr++ = id2;
    *LEstsptr++ = i+1;
    for ( j = 0; j < nb; j++ ) {
      *LEstsptr++ = (logS->LEsts)[i][j];
    }
    *LEstsptr++ = (logS->LIter )[i][0];
    *LEstsptr++ = (logS->Lconvg)[i][0];
    *LEstsptr++ = (logS->Lnbeta)[i][0];
    *LEstsptr++ = (logS->LProb )[i][0];
  }
  *LEptr2 = LEstsptr;

  Modelsptr = *Modelptr2;
  for ( i = 0; i < nt; i++ ) { 
    //    don't add these to output
    //    *Modelsptr++ = id;
    //    *Modelsptr++ = id2;
    //    *Modelsptr++ = i+1;
    for ( j = 0; j < nb; j++ ) {
      *Modelsptr++ = (logS->Models)[i][j];
    }
  }
  *Modelptr2 = Modelsptr;
  return 0;
}

// write optimal P and Q smoothing parameters to output
int toOut( double **ptr, int v1, int v2, int v3, int v4, double v5, double v6 ) {
  double *opt;
  opt = *ptr;

  *opt++ = v1;
  *opt++ = v2;
  *opt++ = v3;
  *opt++ = v4;
  *opt++ = v5;
  *opt++ = v6;

  *ptr = opt;
  return 0;
}
