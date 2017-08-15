// ----------------------------------------------------------------------------
// Header to define Pstruct a structure to hold the conditional probability of
// staying on study for discrete longitudinal data
// ----------------------------------------------------------------------------
#ifndef Pstruct_h_
#define Pstruct_h_

typedef struct Pstruct {
  int NT;         // number of time points (actual)
  int N;          // number of obs in Y

  int Type;       // Type of struct 0 means dropout and
                  // 1 means no dropout

  // Training
  int *Na;        // number of obs in a at each time-point
  double ***a;    // a[t] is Na[t] x 3  

  // Target
  int *Nb;        // number of obs in b at each time-point
  double ***b;    // b[t] is Nb[t] x 3

  // prob of being on study
  double **P;     // P[t] is Nb[t] in length
  double **D1;    // First derivative
  double **D2;    // Second derivative

  int **Posb;     // Posb[t] is Nb[t] in length and contains
                  // the position of b[t][0] in unique list

  //  int **rix;      // row index for look up in b

  int *NTa;       // number onstudy in a at t
  int *NTb;       // number onstudy in b at t

} Pstruct;

// Create a Pstruct but leave the population til later when the data becomes available
Pstruct *Pinit0( int N0, int NT, int size, int type );

// Populate an existing Pstruct -- i.e. refill the a and b tables 
int Pinit1( Pstruct *Xptr, double **Y, int N0, int NT, int start_cut, int stop_cut, int type );

// update an existing Pstruct using the smoothing parameter sigma
int updateP( Pstruct *Xptr, double sigma);

// destruct a Pstruct
int Pdestruct( Pstruct *Xptr );

// print a Pstruct (debug only)
int printPstruct( Pstruct *X );

// This builds the a table for a P object
// caller passes space for output table and it uses SEnv.scratch
int mkPaTablex( double **x, int n, int m, int *Tabrows , double **Table);
#endif
