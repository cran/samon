// ----------------------------------------------------------------------------
// Header to define Qstruct a structure to hold the conditional distribution
// of value Y at time t+1 given Y_t and R_t = 1 in a study with discrete
// longitudinal data.
// ----------------------------------------------------------------------------
#ifndef Qstruct_h_
#define Qstruct_h_

typedef struct Qstruct {
  int NT;           // number of time points (actual)
  int N;            // number of obs in Y

  int Type;         // Type of struct 0 means dropout and
                    // 1 means no dropout

  // Training
  int *Na;          // number of obs in a at each time-point
  double ***a;      // a[t] is Na[t] x 3  

  // Target
  int *Nb;          // number of obs in b at each time-point
  double ***b;      // b[t] is Nb[t] x 3

  int *acnt;        // total in a at t+1
  int *bcnt;        // total in b at t+1

  int *Nr;          // number of rows in q 
  int *Nc;          // number of cols in q

  double **Qr;      // row values in q  (unique from a + b[t][0][])
  double **Qc;      // col values in q  (unique from a + b[t][1][])

  int **Posr;       // position of r value in unique table
  int **Posc;       // position of c value in unique table

  double **Qd;      // col distribution 

  //  int **rix;        // row index for look up
  //  int **cix;        // col index for look up

  double ***Q;      // Nr by Nc Q matrix
  double ***D1;     // Nr by Nc first derivative
  double ***D2;     // Nr by Nc second derivative
  double ***CQ;     // Nr by Nc cumulative Q matrix
  double ***DCQ1;   // Nr by Nc cumulative first derivative
  double ***DCQ2;   // Nr by Nc cumulative second derivative

  double ***TQ;     // Tilted version of Q
  double ***H;      // Used in IF calculations
  double ***IFB;    // Used in IF calculations

} Qstruct;

// Create a Qstruct but leave the population til later when the data becomes available
Qstruct *Qinit0( int N0, int NT, int size, int type );

// Populate an existing Qstruct -- i.e. refill the a and b tables 
int Qinit1( Qstruct *Xptr, double **Y, int N0, int NT, int start_cut, int stop_cut, int type );

// update an existing Qstruct using the smoothing parameter sigma
int updateQ( Qstruct *Xptr, double sigma );

// destruct a Qstruct
int Qdestruct( Qstruct *X );

// print a qstruct (debug only)
int printQstruct( Qstruct *X );

int posc( Qstruct *Xptr, int t, double Y );
int posr( Qstruct *Xptr, int t, double Y );
#endif
