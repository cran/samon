// ----------------------------------------------------------------------------
// Gen_fun: generates data under parametric smoothing.
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"
#include "Qstruct.h"
#include "Gen_fun.h"

#define randNN  103
int First_time = 0;
unsigned long long int randx;
unsigned long long int randarr[randNN];

int posPb( Pstruct *X, int t, double Y );
int posQb( Qstruct *X, int t, double Y );

// ------------------------------------------

int Gen_fun(double **y, int n0, int nt, double **sampy, int Nsamp, double sigp, double sigq ) 
{
  int i, j, ip;
  int t;
  Pstruct *Pptr;
  Qstruct *Qptr;

  int itest;
  double xtest;
  int found;
  int *sampr;

  Pptr = Pinit0( n0, nt, SEnv.nunique, 1);  
  Pinit1( Pptr, y, n0, nt, 0, 0, 1);
  updateP( Pptr, sigp );

  Qptr = Qinit0( n0, nt, SEnv.nunique, 1);
  Qinit1( Qptr, y, n0, nt, 0, 0, 1);
  updateQ( Qptr, sigq );

  sampr = SEnv.Rscratch;

  for ( i = 0; i < Nsamp; i++ ) {
    for ( j = 0; j < nt; j++ ) {
      sampy[i][j] = NAN;
      sampr[j] = -1;
    }
    sampr[0] = 1;

    itest = floor(n0*sgen());
    if (itest > (n0-1)) itest = n0-1;
    sampy[i][0] = y[itest][0];

    for ( t = 1; t < nt; t++ ) {
      if ( sampr[t-1] == 1 ) {
	xtest = sgen();

	ip = posPb( Pptr, t-1, sampy[i][t-1] );
	if ( xtest <= ((Pptr->P)[t-1])[ip] ) sampr[t] = 1;
	else sampr[t] = 0;

	if ( sampr[t] == 1 ) {
	  xtest = sgen();
	  found = 0;
	  j     = 0;
	  ip    = posQb( Qptr, t-1, sampy[i][t-1] );
	  do {
	    if ( xtest <= (Qptr->CQ[t-1])[ ip ][j++]  ) found = 1;
	  } while ( found == 0 && j < (Qptr->Nc[t-1]) );

	  sampy[i][t] = (Qptr->Qc[t-1])[j-1];
	}
      }
      else sampr[t] = 0;
    }
  }
 
  Pdestruct( Pptr );
  Qdestruct( Qptr );

  return 0;
} 

int posPb( Pstruct *X, int t, double Y )
{
  int i;
  if ( Y < ((X->b)[t])[0][0] ) return -1;
  for ( i = 0; i < (X->Nb)[t]; i++ ) {
    if ( Y == ((X->b)[t])[i][0] ) return i;
  }
  return -1;
}

int posQb( Qstruct *X, int t, double Y )
{
  int i;
  if ( Y < ((X->Qr)[t])[0] ) return -1;
  for ( i = 0; i < (X->Nr)[t]; i++ ) {
    if ( Y == ((X->Qr)[t])[i] ) return i;
  }
  return -1;
}

// -------------------------------------------------------------------------------------------

int  seed_sgen(unsigned long long int inseed)
{
  int i;
  static unsigned long long int m = 9223372036854775808llu;
  static unsigned long long int c = 999llu;
  static unsigned long long int a = 892368405llu;

  randx = inseed;
  for ( i = 0; i < randNN; i++ ) {
    randx = ( a * randx + c ) % m;
    randarr[i] = randx;
  }
  return 0;
}

double sgen(void)
{
  int pos;
  double x;
  static unsigned long long int m = 9223372036854775808llu;
  static unsigned long long int c = 999llu;
  static unsigned long long int a = 892368405llu;

  randx = ( a * randx + c ) % m;
  pos  = randx % randNN;
  x = (double) randarr[pos] / (double) m;
  randarr[pos] = randx;
  return x;
}
