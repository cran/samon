// ----------------------------------------------------------------------------
// samon.h header file for samon. 
// ----------------------------------------------------------------------------

#ifndef _samon_h
#define _samon_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "Pstruct.h" 
#include "Qstruct.h" 
#include "Qscratch.h" 
#include "Pscratch.h" 

// -------------------------------------------------------------
   
typedef struct samonEnv {
  int NT;
  int N0;
  int MaxValue;

  double *uvalues;         // unique set of input values sorted low to high
  int nunique;             // Number of values in uvalues

  int seed0;

  double startp;
  double HSigp;
  double startq;
  double HSigq;

  double SmallV;
  double FAconvg;
  double FRconvg;
  double SAconvg;

  double minMainPQ[2];     // main optimal sigmap and sigmaq
  double **minSampPQ;      // sample optimal sigmap and sigmaq;

  int MaxIter;

  double LAlpha;
  double HAlpha;
  double IncAlpha;

  double lb;
  double ub;

  double zeta1;
  double zeta2;

  int Nsamp;
  int NSamples;

  int NParts;
  int **Part;

  Pstruct **Pptrs;
  Qstruct **Qptrs;

  int Nscr;           // dim of scratch items for freeing later

  Pscratch *Pscrsch;
  Qscratch *Qscrsch;

  int *Rscratch;

} samonEnv;

samonEnv SEnv;
#endif



