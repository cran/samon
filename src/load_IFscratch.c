// ----------------------------------------------------------------------------
// load_IF: Load the IF scratch area for use by IF functions
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pscratch.h"
#include "Qscratch.h"
#include "IFscratch.h"
#include "beta_cdf.h"

int load_IF(int N0, int NT, int nuvalues, double **uvalues, int lenTfunMat, double **TfunMat )
{
  int i;
  int MXV, MXT;
// ------------------------------------------------------------

  IFscr.N0 = N0;
  IFscr.NT = NT;

// ------------------------------------------------------------

   IFscr.nuvalues = nuvalues;
   IFscr.uvalues  = uvalues;

   MXV = nuvalues;
   MXT = SEnv.NT + 1;

   IFscr.Pptr = Pinit0( N0, NT, MXV, 1);
   IFscr.Qptr = Qinit0( N0, NT, MXV, 1);

   IFscr.tPMat   = mkMatdz( MXV, MXT           );
   IFscr.tdv     = mkMatdz( MXV, MXT           );
   IFscr.tInt    = mkMatdz( MXV, MXT           );
   IFscr.tBInt   = mkMatdz( MXV, MXT           );
   IFscr.tCInt   = mkMatdz( MXV, MXT           );
   IFscr.tPre    = mkMatdz( MXV, MXT           );
   IFscr.tT2     = mkMatdz( MXV, MXT           );
   IFscr.tB3_3   = mkMatdz( MXV, MXT           );
   IFscr.tCM     = mkMatdz( MXV, MXT           );
   IFscr.Top     = mkMatdz( MXV, MXT           );
   IFscr.Bot     = mkMatdz( MXV, MXT           );

   IFscr.tV      = mkMatiz( MXV, MXT           );

   IFscr.tQ0     = malloc( MXV * sizeof(double));
   IFscr.tEalpha = malloc( MXV * sizeof(double));
   IFscr.trfun   = malloc( MXV * sizeof(double));

   for ( i = 0; i < MXV; i++ ) {
     IFscr.tQ0[i]     = 0.0;
     IFscr.tEalpha[i] = 0.0;
     if ( lenTfunMat > 0 && i < lenTfunMat ) {
       IFscr.trfun[i] = TfunMat[i][1];
     }
     else {
       IFscr.trfun[i] = beta_cdf(( uvalues[i][0] - SEnv.lb) / ( SEnv.ub - SEnv.lb ),SEnv.zeta1,SEnv.zeta2);
     }
   }

   return 0;
}

int free_IF()
{
// ------------------------------------------------------------

  Pdestruct( IFscr.Pptr );
  Qdestruct( IFscr.Qptr );

  freeMatd( IFscr.tPMat );
  freeMatd( IFscr.tdv   );
  freeMatd( IFscr.tInt  );
  freeMatd( IFscr.tBInt );
  freeMatd( IFscr.tCInt );
  freeMatd( IFscr.tPre  );
  freeMatd( IFscr.tT2   );
  freeMatd( IFscr.tB3_3 );
  freeMatd( IFscr.tCM   );
  freeMatd( IFscr.Top   );
  freeMatd( IFscr.Bot   );
  freeMati( IFscr.tV    );

  free( IFscr.tQ0       );
  free( IFscr.tEalpha   );
  free( IFscr.trfun     );

  return 0;
} 
