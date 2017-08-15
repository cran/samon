// ----------------------------------------------------------------------------
// Header to define IFscratch structure.  
// This holds matrices for computing various parts of IF.
// ----------------------------------------------------------------------------

#ifndef IFscratch_h_
#define IFscratch_h_

#include "samon_env.h"
#include "Pstruct.h"
#include "Qstruct.h"

typedef struct IFFscratch {

  double **tPMat;
  //  double **tQMat;
  //  double **tTMat;
  double **tdv;

  double **tInt;
  double **tBInt;
  double **tCInt; 

  double **tPre;
  //  double **tHMat;

  //  double **tBMat;

  double **tT2;
  double **tB3_3;
  double **tCM;

  int MXV;
  int MXT;

  int N0, NT;

  double UU;

  int nuvalues;
  double **uvalues;

  double *tEalpha;
  double *trfun;

  int **tV;
  int **tYY0;

  double *tQ0;
  double **Top;
  double **Bot;

  Pstruct *Pptr;
  Qstruct *Qptr; 

} IFscratch;

IFscratch IFscr;

// load the IFscratch area IFscr defined here
//int load_IF(int N0, int NT, int nuvalues, double **uvalues );
int load_IF(int N0, int NT, int nuvalues, double **uvalues, int lenTfunMat, double **TfunMat );
// and free it
int free_IF();

#endif
