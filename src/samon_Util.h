// ----------------------------------------------------------------------------
// samon_Util.h header file for samon. 
// ----------------------------------------------------------------------------
#ifndef _samon_Util_h
#define _samon_Util_h

#include "samon_logit.h"

int setBase(int **model, int nt, int nb );
int logStoOut(LogisticS *logS,double **LEptr2, int **Modelptr2, int nt, int nb, int id, int id2 );
int toOut( double **ptr, int v1, int v2, int v3, int v4, double v5, double v6 );

#endif
