// ----------------------------------------------------------------------------
// load_env: Load the environment
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pscratch.h"
#include "Qscratch.h"
#include "load_env.h"

// --------------------------------------------------------------------------

int load_env ( 
   int            N0,   int           NT,  
   int         seed0,   double    startp,   double     HSigp,
   double     startq,   double     HSigq,
   double         lb,   double        ub,   double     zeta1,   double     zeta2,
   int        NParts,
   int      NSamples,
   int       MaxIter,   double   FAconvg,   double   FRconvg,   double  SAconvg,
   int       nunique,   double    **uvals,  double    smallv 
	       )
{
  int i;

  int **Part;
  int Nscr;
  Pstruct **Pstarr;
  Qstruct **Qstarr;

  SEnv.NT         = NT;
  SEnv.N0         = N0;
  SEnv.MaxValue   = nunique;
  SEnv.seed0      = seed0;

  SEnv.startp     = startp;
  SEnv.HSigp      = HSigp;

  SEnv.startq     = startq;
  SEnv.HSigq      = HSigq;

  SEnv.MaxIter    = MaxIter;

  if ( smallv > 0.0 ) SEnv.SmallV = smallv;
  else  SEnv.SmallV = 1.0E-2;

  SEnv.FAconvg    = FAconvg;
  SEnv.FRconvg    = FRconvg;
  SEnv.SAconvg    = SAconvg;

  SEnv.lb         = lb;
  SEnv.ub         = ub;
  SEnv.zeta1      = zeta1;
  SEnv.zeta2      = zeta2;

  SEnv.NParts     = NParts;
  SEnv.Nsamp      = N0;
  SEnv.NSamples   = NSamples;

  // --------------------------------------------------------------

  SEnv.nunique    = nunique;
  SEnv.uvalues    = malloc(nunique * sizeof(double) );
  for ( i = 0; i < nunique; i++ ) SEnv.uvalues[i] = uvals[i][0];

  SEnv.minSampPQ  = mkMatd(NSamples,2);

  // --------------------------------------------------------------

  // takes N0, NParts, etc, from the Env 
  partition( &Part, NParts, N0 );
  SEnv.Part       = Part;

  Pstarr = malloc( NParts * sizeof( Pstruct* ));
  SEnv.Pptrs = Pstarr;

  Qstarr = malloc( NParts * sizeof( Qstruct* ));
  SEnv.Qptrs = Qstarr;

  // put P and Q init 0 here
  for ( i = 0; i < NParts; i++ ) {
    SEnv.Pptrs[i] = Pinit0( N0, NT, nunique, 0);
    SEnv.Qptrs[i] = Qinit0( N0, NT, nunique, 0);
  }

  // --------------------------------------------------------------

  // scratch area for sorting etc
  Nscr      = N0 + 1; 
  SEnv.Nscr = Nscr;

  Pscratch *Pscr;
  Pscr = malloc( sizeof( Pscratch ) );

  SEnv.Pscrsch = Pscr;
  (SEnv.Pscrsch)->zmata   = mkMatd( Nscr, 2 );
  (SEnv.Pscrsch)->zmatb   = mkMatd( Nscr, 2 );
  (SEnv.Pscrsch)->utable  = mkMatd( Nscr, 3 );

  Qscratch *Qscr;
  Qscr = malloc( sizeof( Qscratch ) );

  SEnv.Qscrsch = Qscr;
  (SEnv.Qscrsch)->zmata   = mkMatd( Nscr, 3 );
  (SEnv.Qscrsch)->zmatb   = mkMatd( Nscr, 3 );
  (SEnv.Qscrsch)->zmatp0  = malloc( Nscr * sizeof(double) );
  (SEnv.Qscrsch)->zmatp1  = malloc( Nscr * sizeof(double) );
  (SEnv.Qscrsch)->tablep0 = mkMatd( Nscr, 2 );
  (SEnv.Qscrsch)->tablep1 = mkMatd( Nscr, 2 );

  SEnv.Rscratch = malloc( NT * sizeof(int) );
  // --------------------------------------------------------------

  return 0;
}

// Reduce N0 in the environment 
int reload_env ( int N0, double startp, double startq )
{
  int i;
  int NT;
  int NParts;
  int nunique;

  // --------------------------------------------------------------

  int **Part;


  if ( startp < SEnv.HSigp  && SEnv.SmallV < startp ) SEnv.startp = startp;
  if ( startq < SEnv.HSigq  && SEnv.SmallV < startq ) SEnv.startq = startq;

  // --------------------------------------------------------------

  // free and make anew
  NT      = SEnv.NT;
  Part    = SEnv.Part;
  NParts  = SEnv.NParts;
  nunique = SEnv.nunique;

  freeMati( Part );

  partition( &Part, NParts, N0 );
  SEnv.Part       = Part;

  // put P and Q init 0 here
  for ( i = 0; i < NParts; i++ ) {
    Pdestruct( SEnv.Pptrs[i] );
    Qdestruct( SEnv.Qptrs[i] );

    SEnv.Pptrs[i] = Pinit0( N0, NT, nunique, 0);
    SEnv.Qptrs[i] = Qinit0( N0, NT, nunique, 0);
  }
  SEnv.N0         = N0;

  // --------------------------------------------------------------
  return 0;
}


// Free the environment
int free_env ( )
{
  int i;
  int NParts;
  int **Part;

  // --------------------------------------------------------------

  Part     = SEnv.Part;
  NParts   = SEnv.NParts;

  // --------------------------------------------------------------

  if ( SEnv.minSampPQ != NULL ) freeMatd(SEnv.minSampPQ);
  if ( Part           != NULL ) freeMati( Part );

  for ( i = 0; i < NParts; i++ ) {
    Pdestruct( SEnv.Pptrs[i] );
    Qdestruct( SEnv.Qptrs[i] );
  }

  if ( SEnv.Pptrs != NULL ) free(SEnv.Pptrs);
  if ( SEnv.Qptrs != NULL ) free(SEnv.Qptrs);

  // --------------------------------------------------------------

  if (  (SEnv.Pscrsch)->zmata  != NULL )  freeMatd( (SEnv.Pscrsch)->zmata   );
  if (  (SEnv.Pscrsch)->zmatb  != NULL )  freeMatd( (SEnv.Pscrsch)->zmatb   );
  if (  (SEnv.Pscrsch)->utable != NULL )  freeMatd( (SEnv.Pscrsch)->utable  );
  if (   SEnv.Pscrsch          != NULL )  free(      SEnv.Pscrsch           );

  if (  (SEnv.Qscrsch)->zmata  != NULL )  freeMatd( (SEnv.Qscrsch)->zmata   );
  if (  (SEnv.Qscrsch)->zmatb  != NULL )  freeMatd( (SEnv.Qscrsch)->zmatb   );
  if (  (SEnv.Qscrsch)->zmatp0 != NULL )  free(     (SEnv.Qscrsch)->zmatp0  );
  if (  (SEnv.Qscrsch)->zmatp1 != NULL )  free(     (SEnv.Qscrsch)->zmatp1  );
  if (  (SEnv.Qscrsch)->tablep0!= NULL )  freeMatd( (SEnv.Qscrsch)->tablep0 );
  if (  (SEnv.Qscrsch)->tablep1!= NULL )  freeMatd( (SEnv.Qscrsch)->tablep1 );
  if (   SEnv.Qscrsch          != NULL )  free(      SEnv.Qscrsch           );

  if (   SEnv.Rscratch         != NULL )  free(      SEnv.Rscratch          );

  if (   SEnv.uvalues          != NULL )  free(      SEnv.uvalues           );
  return 0;
}
