// ----------------------------------------------------------------------------
// header defining a Pscratch structure.  
// Such structures are embedded in the environment and serve in the creation
// of Pstructs.
// ----------------------------------------------------------------------------
#ifndef Pscratch_h_
#define Pscratch_h_
typedef struct Pscratch {
  double **zmata;       // Nscr = (N0+1) by 2 
  double **zmatb;       // Nscr = (N0+1) by 2

  double **utable;      // Nscr = (N0+1) by 3
} Pscratch;
#endif
