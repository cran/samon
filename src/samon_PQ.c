// ----------------------------------------------------------------------------
// samon_PQ: Entry point for samon in R, generate a bootstrap sample 
// ----------------------------------------------------------------------------
// 0. Initialize environment
// 1. Popt, Qopt for input data 
// 3. Return the P and Q matrices 
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_IF ( 
//  int               *N0,   Number of observations
//  int               *NT,   Number of time-points
//  double           *Mat,   *N0 by *NT input matrix
//
//  double      *outPmatM,   Output matrix for main P smoothing parameter
//  double      *outQmatM,   Output matrix for main Q smoothing parameter
//
//  int            *seed0,   initial seed 
//  double        *startp,   Start value to use for P smoothing parameter
//  double         *HSigp,   High value for P smoothing
//  double        *startq,   Start value tor use for Q smoothing parameter
//  double         *HSigq,   High value for Q smoothing
//
//  int           *NParts,   Number of partitions
//  int         *NSamples,   Dummies
//
//  int          *MaxIter,   Convergence criteria for optimization
//  double       *FAconvg,
//  double       *FRconvg,
//  double       *SAconvg,
//
//  double          *Pmat,   output P matrix
//  double          *Qmat,   output Q matrix
//  int          *lenTfun,   Tilt function nrow
//  double          *Tfun    Tilt matrix (*lenTfun by 2)
//             )
// ------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"
#include "Qstruct.h"
#include "IFscratch.h"
#include "samon_opt.h"
#include "load_env.h"
#include "Gen_fun.h"
#include "IF_fun_beta.h"
#include "samon_Util.h"

#include "unique.h"
// ------------------------------------------------------------------------
void samon_PQ ( 
    int            *N0,    int           *NT,   double       *Mat, 
    double   *outPmatM,   double  *outQmatM,

    int         *seed0,
    double     *startp,   double     *HSigp,
    double     *startq,   double     *HSigq,

    int        *NParts,
    int      *NSamples,

    int       *MaxIter,
    double    *FAconvg,   double   *FRconvg,    double  *SAconvg,
    double       *Pmat,   double      *Qmat,
    int       *lenTfun,   double      *Tfun
		)
{
  int i, j, t;
  double *ptr;
  int rc;

  double *optr;
  double *optPmatp, *optQmatp;

  double **TfunMat = NULL;

  // main 
  int    MiterP, MiterQ, *MiterPp, *MiterQp;
  double MxminP, MxminQ, *MxminPp, *MxminQp;
  double MfminP, MfminQ, *MfminPp, *MfminQp;

  Pstruct *Pptr;
  Qstruct *Qptr; 

  // --------------------------------------------------------------
  double **dMat;
  double MaximumValue, MinimumValue;
  // --------------------------------------------------------------

  dMat    = mkMatd(*N0, *NT);

  ptr = Mat;
  if ( isnan(*ptr) ) return;

  MaximumValue = *ptr;
  MinimumValue = *ptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      dMat[i][j] = *ptr++;
      if ( !isnan(dMat[i][j]) && dMat[i][j] > MaximumValue ) MaximumValue = dMat[i][j];
      if ( !isnan(dMat[i][j]) && dMat[i][j] < MinimumValue ) MinimumValue = dMat[i][j];
    }
  }

  double **uvals;
  double smallinc;
  int count;
  uvals = uniqueVal( dMat, *N0, *NT, &count );

  smallinc = 0.1;
  for ( i = 1; i < count; i++ ) {
    if (i == 1) smallinc = uvals[i][0] - uvals[i-1][0];
    else {
      if ((uvals[i][0] - uvals[i-1][0]) < smallinc) smallinc = uvals[i][0] - uvals[i-1][0];
    }
  }

  if ( *lenTfun > 0 ) {
    TfunMat = mkMatd(*lenTfun,2);
    ptr = Tfun;
    for ( j = 0; j < 2; j++ ) {
      for ( i = 0; i < *lenTfun; i++ ) {
	//    for ( j = 0; j < 2; j++ ) {
	//      for ( i = 0; i < *lenTfun; i++ ) {
	TfunMat[i][j] = *ptr++;
      }
    }
  }

  rc = load_env ( *N0, *NT, *seed0, *startp, *HSigp, *startq, *HSigq,
		  0, 0, 1, 1,
                  *NParts,
                  *NSamples,
                  *MaxIter, *FAconvg, *FRconvg, *SAconvg, 
                  count, uvals, smallinc/4.0);

  rc = load_IF(*N0,*NT,count, uvals, *lenTfun, TfunMat );

  // 1.0
  // --------------------------------------------------------------
  MiterPp  = &MiterP;
  MxminPp  = &MxminP;
  MfminPp  = &MfminP; 
  MiterQp  = &MiterQ;
  MxminQp  = &MxminQ;
  MfminQp  = &MfminQ; 

  optPmatp = outPmatM;
  optQmatp = outQmatM;

  // P opt
  rc =  Popt(dMat, MiterPp, MxminPp, MfminPp);
  rc =  toOut( &optPmatp, 0, 0, rc, MiterP, MxminP, MfminP);

  // Q opt
  rc =  Qopt(dMat, MiterQp, MxminQp, MfminQp);
  rc =  toOut( &optQmatp, 0, 0, rc, MiterQ, MxminQ, MfminQ);

  Pptr = IFscr.Pptr;
  rc   = Pinit1(Pptr, dMat, *N0, *NT, 0, 0, 1);
  rc   = updateP( Pptr, MxminP );

  Qptr = IFscr.Qptr;
  rc   = Qinit1(Qptr, dMat, *N0, *NT, 0, 0, 1);
  rc   = updateQ( Qptr, MxminQ );

  optr = Pmat;
  for ( t = 0; t < Pptr->NT; t++ ) {
    for ( i = 0; i < (Pptr->Nb)[t]; i++ ) {
      *optr++  = t;
      *optr++  = i;
      *optr++  =  (Pptr->b[t])[i][0];
      *optr++  =  (Pptr->b[t])[i][1];
      *optr++  =  (Pptr->b[t])[i][2];
      *optr++  =  (Pptr->P[t])[i];
    }
  }

  optr = Qmat;
  for ( t = 0; t < Qptr->NT; t++ ) {
    for ( i = 0; i < (Qptr->Nr)[t]; i++ ) {
      for ( j = 0; j < (Qptr->Nc)[t]; j++ ) {
	*optr++  = t;
	*optr++  = i;
	*optr++  = j;
	*optr++  =  (Qptr->Qr[t])[i];
	*optr++  =  (Qptr->Qc[t])[j];
	*optr++  =  (Qptr->Q[t])[i][j];
      }
    }
  }
  // --------------------------------------------------------------
  freeMatd(dMat    );
  rc = free_IF();
  rc = free_env();
  freeMatd(uvals);
  if ( *lenTfun > 0 ) freeMatd( TfunMat );
  return;
}
