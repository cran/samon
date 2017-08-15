// ----------------------------------------------------------------------------
// Optimize the P smoothing parameter sigma_p via the loss function lossP
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"
#include "samon_opt.h"

int Popt(double **y, int *iter, double *optx, double *optfn) {
  int i;
  int N0, NT;
  int minrc;
  int NParts;

  N0     = SEnv.N0;
  NT     = SEnv.NT;
  NParts = SEnv.NParts;

  for ( i = 0; i < NParts; i++ ) {
    Pinit1( SEnv.Pptrs[i], y, N0, NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);
  }

  minrc = Pmin(NParts, iter, optx, optfn);

  return minrc; 
}

// ---------------------------------------------------------------------------

int Pmin(int NParts, int *iter, double *optx, double *optfn) {
  int rc, i;
  double x0, stp;
  double der1, der2;
  double *der1ptr, *der2ptr;
  double xh;
  double der1h, der2h;
  int MaxIter;
  double BigValue, SmallValue;
  double StartValue;
  double distance;

  double fn0, fn, fnB, fn1, fnH;
  double FAconvg, FRconvg, SAconvg;
  double f[7], x[7], d1[7], d2[7];

  MaxIter    = SEnv.MaxIter;
  FAconvg    = SEnv.FAconvg;
  FRconvg    = SEnv.FRconvg;
  SAconvg    = SEnv.SAconvg;

  SmallValue = SEnv.SmallV;
  BigValue   = SEnv.HSigp;
  StartValue = SEnv.startp;

  der1 = 0.0;
  der2 = 0.0;

  der1ptr = &der1;
  der2ptr = &der2;

  der1h = 0.0;
  der2h = 0.0;

    if ( StartValue < SmallValue ) StartValue = 0.7 * SmallValue + 0.3 * BigValue ;  
    if ( StartValue >   BigValue ) StartValue = 0.5 * SmallValue + 0.5 * BigValue ;
    distance = StartValue - SmallValue;
    if ((BigValue - StartValue) < distance ) distance = BigValue - StartValue;

  x0 = StartValue;
  if ( x0 < SmallValue ) x0 = SmallValue;
  fn0 = lossP(x0, NParts, der1ptr, der2ptr);

  int mark;
  mark  = 0;
  xh    = x0;
  fnH   = fn0;
  der1h = der1;
  der2h = der2;
  for ( i = 0; i < 6; i++ ) {
    x[i] = StartValue - (3.5-i) * distance / 10.0;
    f[i] = lossP(x[i], NParts, &(d1[i]), &(d2[i]));
    if ( f[i] < fnH && (fabs(d2[i]) > 1.0E-50) ) {
      mark  = 1;
      xh    = x[i];
      fnH   = f[i];
      der1h = d1[i];
      der2h = d2[i];
    }
  }

  if ( mark == 1) {
    x0   = xh;
    fn0  = fnH;
    der1 = der1h;
    der2 = der2h;
  }

  *iter = 0;
  fn1   = fn0;
  rc    = 0;
  stp   = 0.5;
  fn    = fn0;
  do {
    if ( fabs(der2) < 1.0E-50 ) {
      // derivative too small to take step 
      rc = 3;
    }
    else {
      stp = ( der1 / der2 );

      if ( (x0-stp) < SmallValue ) {
	do {
	  stp = stp / 10;
	} while ( (x0-stp) < SmallValue );
	stp = stp / 2;
	x0 = x0 - stp;
	fn = lossP(x0, NParts, der1ptr, der2ptr);
      }
      else if ( (x0 - stp) > BigValue ) {
	// stop if you go beyond BigValue
	x0 = BigValue;
	fn = lossP(x0, NParts, der1ptr, der2ptr);
	rc = 5;
      }
      else {
	x0 = x0 - stp;
	fn = lossP(x0, NParts, der1ptr, der2ptr);
      }
    if ( fabs(  fn - fn1 ) < FAconvg ) rc = 1;
    if ( fabs( (fn - fn1) / fabs(fn+fn1) ) < FRconvg ) rc = 2;
    }
    fn1 = fn;
  } while ( (*iter)++ < MaxIter && fabs(stp) > SAconvg && rc == 0);  

  if ( rc == 0 && *iter >= MaxIter ) rc = 4;

  fnB = lossP(BigValue,   NParts, der1ptr, der2ptr);

  if ( fnB < fn ) {
    x0 = BigValue;
    fn = fnB;
    rc = 7;
  } 

  *optx  = x0;
  *optfn = fn; 
  if ( *iter == 0 ) {
    x0 = StartValue;
    *optfn = fn0;
  }

  return rc;
}

double lossP( double sigma, int NParts, double *deriv1, double *deriv2 )
{
  int p, t, i;
  int nt, nb;
  double n1, n2;
  double pv, lossp, lossi;
  double der1, der2;
  double d1, d2;
  double nCut;

  nt = (SEnv.Pptrs[0])->NT;

  lossp = 0.0;
  der1  = 0.0;
  der2  = 0.0;
  for ( p = 0; p < NParts; p++ ) {
    updateP(SEnv.Pptrs[p], sigma );
    nCut = SEnv.Part[p][1] - SEnv.Part[p][0] + 1;
    lossi = 0;
    for ( t = 0; t < (nt-1) ; t++ ) {
      nb = (SEnv.Pptrs[p]->Nb)[t];

      for( i = 0; i < nb; i++ ) {
	pv = ((SEnv.Pptrs[p]->P )[t])[i];
        d1 = ((SEnv.Pptrs[p]->D1)[t])[i];
        d2 = ((SEnv.Pptrs[p]->D2)[t])[i];
        n1 = ((SEnv.Pptrs[p]->b)[t])[i][1];
        n2 = ((SEnv.Pptrs[p]->b)[t])[i][2];

	lossi = lossi + ( n2/nCut ) * ( 1 - pv ) * ( 1 - pv );
	lossi = lossi + ( n1/nCut ) * ( 0 - pv ) * ( 0 - pv );

	der1  = der1  + ( n2/nCut ) * ( -2.0 ) * ( 1 - pv ) * d1
	              + ( n1/nCut ) * (  2.0 ) * (     pv ) * d1;
        der2  = der2  + ( n2/nCut ) * (  2.0 ) * ( pow(d1,2) - ( 1 - pv) * d2 )
	              + ( n1/nCut ) * (  2.0 ) * ( pow(d1,2) + (     pv) * d2 );

      }
    }
    lossp = lossp + lossi;
  }
  *deriv1 = der1;
  *deriv2 = der2;

  return lossp;
}
