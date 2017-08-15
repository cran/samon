// ----------------------------------------------------------------------------
// Optimize the Q smoothing parameter sigma_q via the loss function lossQ
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Qstruct.h"
#include "samon_opt.h"

int Qopt(double **y, int *iter, double *optx, double *optfn) {
  int i;
  int minrc;
  int N0, NT; 
  int NParts;

  N0     = SEnv.N0;
  NT     = SEnv.NT;
  NParts = SEnv.NParts;

  for ( i = 0; i < NParts; i++ ) {
    Qinit1( SEnv.Qptrs[i], y, N0, NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);
  }

  minrc = Qmin(NParts, iter, optx, optfn);

  return minrc; 
}

// ---------------------------------------------------------------------------

int Qmin(int NParts, int *iter, double *optx, double *optfn) {
  int rc;
  double x0, stp;
  double der1, der2;
  double *der1ptr, *der2ptr;

  int MaxIter;
  double BigValue, SmallValue;
  double StartValue;

  double fn0, fn, fnB, fn1, fnS;
  double FAconvg, SAconvg;

  MaxIter    = SEnv.MaxIter;
  FAconvg    = SEnv.FAconvg;
  SAconvg    = SEnv.SAconvg;

  SmallValue = SEnv.SmallV;
  BigValue   = SEnv.HSigq;
  StartValue = SEnv.startq;

  der1 = 0.0;
  der2 = 0.0;

  der1ptr = &der1;
  der2ptr = &der2;

  x0 = StartValue;
  if ( x0 < SmallValue ) x0 = SmallValue;
  fn0 = lossQ(x0, NParts, der1ptr, der2ptr);

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
      if ( fabs(stp) > 2*SmallValue ) {
	if ( stp < 0 ) stp = -2 * SmallValue;
	else stp = 2 * SmallValue;
      }

      if ( (x0-stp) < SmallValue ) {
	do {
	  stp = stp / 10;
	} while ( (x0-stp) < SmallValue );
	stp = stp / 2;
	x0 = x0 - stp;
	fn = lossQ(x0, NParts, der1ptr, der2ptr);
      }
      else if ( (x0 - stp) > BigValue ) {
	// stop if you go beyond BigValue
	x0 = BigValue;
	fn = lossQ(x0, NParts, der1ptr, der2ptr);
	rc = 5;
      }
      else {
	x0 = x0 - stp;
	fn = lossQ(x0, NParts, der1ptr, der2ptr);
      }
    if ( fabs( fn - fn1 ) < FAconvg ) rc = 1;
    if ( fabs( (fn - fn1) / fabs(fn+fn1) ) < FAconvg ) rc = 2;
    }
    fn1 = fn;
  } while ( (*iter)++ < MaxIter && fabs(stp) > SAconvg && rc == 0);  

  if ( rc == 0 && *iter >= MaxIter ) rc = 4;

  fnS = lossQ(SmallValue, NParts, der1ptr, der2ptr);
  fnB = lossQ(BigValue,   NParts, der1ptr, der2ptr);
  if ( fnS < fn ) {
    x0 = SmallValue;
    fn = fnS;
    rc = 6;
  } 

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

// ---------------------------------------------------------------------------

double lossQ(double sigma, int NParts, double *dlossqptr, double *d2lossqptr)
{
  int i, t, j, p;
  int nt, nb, nc;
  int pos;
  double thisY0, thisY, nn;
  double lossq;
  double dlossq, d2lossq;
  double losst;
  double dlosst, d2losst;
  double nCut;
  
  nt = (SEnv.Qptrs[0])->NT;

  lossq   = 0.0;
  dlossq  = 0.0;
  d2lossq = 0.0;
  for ( p = 0; p < NParts; p++ ) {

    updateQ(SEnv.Qptrs[p], sigma);
    nCut = SEnv.Part[p][1] - SEnv.Part[p][0] + 1;

    for ( t=0; t < (nt-1); t++ ) {

      nb = (SEnv.Qptrs[p])->Nb[t];
      nc = (SEnv.Qptrs[p])->Nc[t];

      losst   = 0.0;
      dlosst  = 0.0;
      d2losst = 0.0;  
      for ( i = 0; i < nb; i++ ) {
	thisY0 = ((SEnv.Qptrs[p]->b)[t])[i][0];
	thisY  = ((SEnv.Qptrs[p]->b)[t])[i][1];
	nn     = ((SEnv.Qptrs[p]->b)[t])[i][2];
	//	pos    = ((SEnv.Qptrs[p]->rix)[t])[ thisY0 ];
	pos   = posr(SEnv.Qptrs[p], t, thisY0 );

	for ( j = 0; j < nc; j++ ) { 
	  if ( thisY <= (((SEnv.Qptrs[p]->Qc)[t]))[j] ) {
	    losst   = losst   + nn * ((SEnv.Qptrs[p]->Qd)[t])[j]       *   pow( ( 1.0 - ((SEnv.Qptrs[p]->CQ  )[t])[pos][j] ), 2.0 );
	    dlosst  = dlosst  + nn * ((SEnv.Qptrs[p]->Qd)[t])[j] * 2.0 *        ( 1.0 - ((SEnv.Qptrs[p]->CQ  )[t])[pos][j] )*(  - ((SEnv.Qptrs[p]->DCQ1)[t])[pos][j]); 
	    d2losst = d2losst + nn * ((SEnv.Qptrs[p]->Qd)[t])[j] * 2.0 * ( pow( (     - ((SEnv.Qptrs[p]->DCQ1)[t])[pos][j] ), 2.0 ) + ( 1.0 - ((SEnv.Qptrs[p]->CQ)[t])[pos][j] )*(  - ((SEnv.Qptrs[p]->DCQ2)[t])[pos][j])); 
	  }
	  else {
	    losst   = losst   + nn * ((SEnv.Qptrs[p]->Qd)[t])[j] *         pow( ((SEnv.Qptrs[p]->CQ  )[t])[pos][j], 2.0 );
	    dlosst  = dlosst  + nn * ((SEnv.Qptrs[p]->Qd)[t])[j] * 2.0 *      ( ((SEnv.Qptrs[p]->CQ  )[t])[pos][j]      )*( ((SEnv.Qptrs[p]->DCQ1)[t])[pos][j]) ; 
	    d2losst = d2losst + nn * ((SEnv.Qptrs[p]->Qd)[t])[j] * 2.0 * ( pow( ((SEnv.Qptrs[p]->DCQ1)[t])[pos][j], 2.0 ) + ( ((SEnv.Qptrs[p]->CQ)[t])[pos][j] )*( ((SEnv.Qptrs[p]->DCQ2)[t])[pos][j])) ; 
	  }
	}

      }
      lossq   = lossq   + losst   / nCut;
      dlossq  = dlossq  + dlosst  / nCut;
      d2lossq = d2lossq + d2losst / nCut;
    }
  }
  *dlossqptr  = dlossq;
  *d2lossqptr = d2lossq;
  return lossq;
}
