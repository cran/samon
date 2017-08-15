// ----------------------------------------------------------------------------
// samon_evalIM: Entry point for samon in R -- intermittent missing data
// ----------------------------------------------------------------------------
// 0. Fill in intermittent missing data
// 1. Popt, and Qopt for input data 
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_evalIM ( 
//  int               *N0,   Number of observations
//  int               *NT,   Number of time-points
//  double           *Mat,   *N0 by *NT input matrix
//
//  int           *nmodel,   size of model (number of betas)
//  int          *inmodel,   input model *NT * *nmodel
//
//  double          *FMat,   Filled in main data
//  double        *LEstsM,   Logistic results for main data
//  int          *ModelsM,   Models for main data
//
//  double        *outmat,   Output matrix 
//
//  int           *nsigma,   number of sigmas
//  double     *sigmaList,   list of sigmas
//
//  int            *seed1,
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
#include "Gen_fun.h"
#include "samon_logit.h"
#include "samon_FillRemove.h"
#include "samon_Util.h"

#include "unique.h"

// -----------------------------------------------------------------------
void samon_evalIM ( 
   int            *N0,  int            *NT,  double        *Mat, 
   int        *nmodel,  int       *inmodel,
   double       *FMat,  double     *LEstsM,  int      *ModelsM,
   double     *outmat,
   int        *nsigma,  double  *SigmaList,
   int         *seed1,
   int        *NParts,
   int          *PQtp
   )
{
  int i, j, nb;
  double *dptr;

  double *outmatp;
  double Result;
  double *LEstsptr;
  int    *Modelsptr;
  int    **BaseModel;
  int    *iptr;

  double MaximumValue, MinimumValue; 

  // --------------------------------------------------------------
  double **iMat;
  double **dMat;

  LogisticS *logS;

  double deriv1, deriv2;
  // --------------------------------------------------------------
  iMat = mkMatd(*N0, *NT);
  dMat = mkMatd(*N0, *NT);

  dptr = Mat;
  if ( isnan(*dptr) ) return;

  MaximumValue = *dptr;
  MinimumValue = *dptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      //for ( i = 0; i < *N0; i++ ) { 
      //    for ( j = 0; j < *NT; j++ ) {
      dMat[i][j] = *dptr++;
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


  // 0.0  Fill in intermittent missing data
  // --------------------------------------------------------------
  seed_sgen(*seed1);

  nb = *nmodel;
  BaseModel = mkMati(*NT,nb);

  iptr = inmodel;
  for ( j = 0; j <  nb; j++ ) {
    for ( i = 0; i < *NT; i++ )
  //  for ( i = 0; i < *NT; i++ )
  //    for ( j = 0; j <  nb; j++ ) {
      BaseModel[i][j] = *iptr++;
  }

  logS  = initLogisticS( dMat, *N0, *NT, nb, 25, 10E-10, 10E-10 );
  cpMati( logS->Models, BaseModel, *NT, nb ); 

  FillIn( nb, logS );

  cpMati( BaseModel, logS->Models, *NT, nb );

  LEstsptr  = LEstsM;
  Modelsptr = ModelsM;
  logStoOut( logS, &LEstsptr, &Modelsptr, *NT, nb, 0, 1 );

  cpMatd(iMat,logS->Data,*N0,*NT); 

  dptr = FMat;
  for ( i = 0; i < *N0; i++ ) { 
    *dptr++ = 0;
    *dptr++ = 1;
    for ( j = 0; j < *NT; j++ ) {
      iMat[i][j] = (logS->Data)[i][j];
      *dptr++ = iMat[i][j];
    }
  }
  distructLogisticS(logS);

  // --------------------------------------------------------------
  if ( *PQtp == 1 ) {
    for ( i = 0; i < *NParts; i++ ) 
      Pinit1( SEnv.Pptrs[i], iMat, *N0, *NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);

    outmatp = outmat;
    for ( i = 0; i < *nsigma; i++ ) {
      Result = lossP( SigmaList[i], *NParts, &deriv1, &deriv2 );
      *outmatp++ = SigmaList[i];
      *outmatp++ = Result;
    }
  } 
  else {
    for ( i = 0; i < *NParts; i++ ) 
      Qinit1( SEnv.Qptrs[i], iMat, *N0, *NT, SEnv.Part[i][0], SEnv.Part[i][1], 0);

    outmatp = outmat;
    for ( i = 0; i < *nsigma; i++ ) {
      Result = lossQ( SigmaList[i], *NParts, &deriv1, &deriv2 );
      *outmatp++ = SigmaList[i];
      *outmatp++ = Result;
    }
  }
// -----------------------------------------------------------------------

  freeMatd( iMat      );
  freeMatd( dMat      );
  freeMati( BaseModel );
  freeMatd( uvals     );
  free_env();
  return;
}
