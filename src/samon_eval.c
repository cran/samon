// ----------------------------------------------------------------------------
// samon_eval: Entry point for samon in R -- intermittent missing data
// ----------------------------------------------------------------------------
// 0. initiate samon
// 1. Calculate Ploss or Qloss for input data 
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_eval ( 
//  int               *N0,   Number of observations
//  int               *NT,   Number of time-points
//  double           *Mat,   *N0 by *NT input matrix
//
//  double        *outmat,   Output matrix 
//
//  int           *nsigma,   number of sigmas
//  double     *sigmaList,   list of sigmas
//
//  int           *NParts,   Number of patritions
//  int             *PQtp    1 means sigma P, 2 means sigma Q 
//             )
// ------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"
#include "Qstruct.h"
#include "samon_opt.h"
#include "load_env.h"

#include "unique.h"
// -----------------------------------------------------------------------
void samon_eval ( 
   int            *N0,    int            *NT,  double        *Mat, 
   double     *outmat,
   int        *nsigma, double    *SigmaList,
   int        *NParts,    int          *PQtp
    )
{
  int i, j;
  double *ptr;

  double *outmatp;
  double Result;

  // --------------------------------------------------------------
  double **dMat;
  double MaximumValue, MinimumValue;

  double deriv1, deriv2;
  // --------------------------------------------------------------

  dMat    = mkMatd(*N0, *NT);

  ptr = Mat;
  if ( isnan(*ptr) ) return;

  MaximumValue = *ptr;
  MinimumValue = *ptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      //for ( i = 0; i < *N0; i++ ) { 
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

  load_env ( *N0, *NT, 1, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0,
                  *NParts, 0,
                  0, 0.0, 0.0, 0.0,
                  count, uvals, smallinc/4.0);

  // --------------------------------------------------------------

  if ( *PQtp == 1 ) {
    for ( i = 0; i < *NParts; i++ ) {
      Pinit1( SEnv.Pptrs[i], dMat, *N0, *NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);
    }

    outmatp = outmat;
    for ( i = 0; i < *nsigma; i++ ) {
      Result = lossP( SigmaList[i], *NParts, &deriv1, &deriv2 );
      *outmatp++ = SigmaList[i];
      *outmatp++ = Result;
    }
  } 
  else {
    for ( i = 0; i < *NParts; i++ ) {
      Qinit1( SEnv.Qptrs[i], dMat, *N0, *NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);
    }

    outmatp = outmat;
    for ( i = 0; i < *nsigma; i++ ) {
      Result = lossQ( SigmaList[i], *NParts, &deriv1, &deriv2 );
      *outmatp++ = SigmaList[i];
      *outmatp++ = Result;
    }
  }

// -----------------------------------------------------------------------

  freeMatd(uvals);
  freeMatd(dMat);
  free_env();
  return;
}
