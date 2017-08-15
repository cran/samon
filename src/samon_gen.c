// ----------------------------------------------------------------------------
// samon_gen: Entry point for samon in R
// ----------------------------------------------------------------------------
// 0. Initialize environment
// 1. Popt, Qopt and IF for input data 
// 2. Generate NSamples ( if *NSamples > 0 )
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_gen ( 
//  int               *N0,   Number of observations
//  int               *NT,   Number of time-points
//  double           *Mat,   *N0 by *NT input matrix
//
//  double      *outPmatM,   Output matrix for main P smoothing parameter
//  double      *outQmatM,   Output matrix for main Q smoothing parameter
//
//  int            *seed0,   Seed for generating bootstraps
//  double        *startp,   Start value to use for P smoothing parameter
//  double         *HSigp,   High value for P smoothing
//  double        *startq,   Start value tor use for Q smoothing parameter
//  double         *HSigq,   High value for Q smoothing
//
//  int           *NParts,   Number of partitions
//  int         *NSamples,   Number of bootstrap samples to create
//
//  int          *MaxIter,   Convergence criteria for optimization
//  double       *FAconvg,
//  double       *FRconvg,
//  double       *SAconvg,
//
//  double        *Sample    The returned sample
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
void samon_gen ( 
  int               *N0,   int               *NT,    double          *Mat, 

  double      *outPmatM,   double      *outQmatM,

  int            *seed0,
  double        *startp,   double         *HSigp,
  double        *startq,   double         *HSigq,
  int           *NParts,
  int         *NSamples,
  int          *MaxIter,
  double       *FAconvg,   double       *FRconvg,   double       *SAconvg,
  double        *Sample
                )
{
  int rc, i, j ;
  double *ptr;

  double *optPmatp, *optQmatp;

  // main 
  int    MiterP, MiterQ, *MiterPp, *MiterQp;
  double MxminP, MxminQ, *MxminPp, *MxminQp;
  double MfminP, MfminQ, *MfminPp, *MfminQp;

  // --------------------------------------------------------------
  double **dMat;
  double **sMat = NULL;
  double MaximumValue, MinimumValue;
  // --------------------------------------------------------------

  dMat    = mkMatd(*N0, *NT);

  ptr = Mat;
  if ( isnan(*ptr) ) return;

  MaximumValue = *ptr;
  MinimumValue = *ptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      //  for ( i = 0; i < *N0; i++ ) { 
      //    for ( j = 0; j < *NT; j++ ) {
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

  rc = load_env ( *N0, *NT, *seed0, *startp, *HSigp, *startq, *HSigq,
                  0, 1, 1, 1,
                  *NParts, *NSamples,
                  *MaxIter, *FAconvg, *FRconvg, *SAconvg, 
                  count, uvals, smallinc/4.0);

  rc = load_IF(*N0,*NT,count, uvals, 0, NULL );

  // 1.0  The main data
  // --------------------------------------------------------------

  MiterPp = &MiterP;
  MxminPp = &MxminP;
  MfminPp = &MfminP; 
  MiterQp = &MiterQ;
  MxminQp = &MxminQ;
  MfminQp = &MfminQ; 

  optPmatp = outPmatM;
  optQmatp = outQmatM;

  // P opt
  rc =  Popt(dMat, MiterPp, MxminPp, MfminPp);
  rc =  toOut( &optPmatp, 0, 0, rc, MiterP, MxminP, MfminP);
  SEnv.minMainPQ[0] = MxminP;

  // Q opt
  rc =  Qopt(dMat, MiterQp, MxminQp, MfminQp);
  rc =  toOut( &optQmatp, 0, 0, rc, MiterQ, MxminQ, MfminQ);
  SEnv.minMainPQ[1] = MxminQ;

  Pstruct *Pptr;
  Qstruct *Qptr; 

  Pptr = IFscr.Pptr;
  rc   = Pinit1(Pptr, dMat, *N0, *NT, 0, 0, 1);
  rc   = updateP( Pptr, MxminP );

  Qptr = IFscr.Qptr;
  rc   = Qinit1(Qptr, dMat, *N0, *NT, 0, 0, 1);
  rc   = updateQ( Qptr, MxminQ );

  // IF
  //  optr     = alphamatM;
  //  ifiptr   = &IFvalsM;
  //  rc = IF_fun(dMat, MxminP, MxminQ, 0, 0, &optr, *Nalpha, alphaList, *retIFiM, ifiptr);

  // 2.0  The samples
  // --------------------------------------------------------------
  if ( *NSamples > 0 ) { 

    // Generate
    seed_sgen(SEnv.seed0);
 
    sMat = mkMatd(*NSamples, *NT);
    rc = Gen_fun( dMat, *N0, *NT, sMat, *NSamples, MxminP , MxminQ );

    ptr = Sample;
    for ( i = 0; i < *NSamples; i++ ) {
      for ( j = 0; j < *NT; j++ ) {
	*ptr++ = sMat[i][j];
      }
    }  

  }

  freeMatd( dMat );

  if ( *NSamples > 0 ) {
    freeMatd(sMat);
  }
  rc = free_IF();
  rc = free_env();
  freeMatd(uvals);
  return;
}
// --------------------------------------------------------------
