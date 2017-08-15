// ----------------------------------------------------------------------------
// header defining a Qscratch structure.  Such structures are embedded in the 
// environment and serve in the creation of Qstructs
// ----------------------------------------------------------------------------
#ifndef Qscratch_h_
#define Qscratch_h_
typedef struct Qscratch {
  double **zmata;       // N0 by 3 
  double **zmatb;       // N0 by 3

  double *zmatp0;       // N0
  double *zmatp1;       // N0

  double **tablep0;     // N0 by 2
  double **tablep1;     // N0 by 2    
} Qscratch;
#endif
