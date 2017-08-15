// ----------------------------------------------------------------------------
// samon_genIM: Entry point for samon in R -- intermittent missing data
// ----------------------------------------------------------------------------
// 0. Initialize environment
// 1. Fill intermittent data n times and for each fill in
// 2. compute Popt, Qopt and IF 
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_ngenIMIF ( 
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
//  double      *outPmatM,   Output matrix for main P smoothing parameter
//  double      *outQmatM,   Output matrix for main Q smoothing parameter
//
//  int            *seed0,   Seed for generating bootstraps (not used here)
//  int            *seed1,   Seed filling in intermittent missing data
//  double        *startp,   Start value to use for P smoothing parameter
//  double         *HSigp,   High value for P smoothing
//  double        *startq,   Start value tor use for Q smoothing parameter
//  double         *HSigq,   High value for Q smoothing
//
//  int           *NParts,   Number of partitions
//  int         *NSamples,   Number of samples to make
//
//  int          *MaxIter,   Convergence criteria for optimization
//  double       *FAconvg,
//  double       *FRconvg,
//  double       *SAconvg,
//
//  int           *NFills,   Number of fillins to do
//
//  double        *Sample    The output sample
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
#include "samon_logit.h"
#include "samon_FillRemove.h"
#include "samon_Util.h"

#include "unique.h"

// ------------------------------------------------------------------------
void samon_genIM ( 
   int            *N0,  int           *NT,  double       *Mat, 
   int        *nmodel,  int      *inmodel,
   double       *FMat,  double    *LEstsM,  int      *ModelsM,
   double   *outPmatM,  double  *outQmatM,
   int         *seed0,  int        *seed1,
   double     *startp,  double     *HSigp,
   double     *startq,  double     *HSigq,
   int        *NParts,
   int      *NSamples, 
   int       *MaxIter,  double   *FAconvg,  double   *FRconvg,  double  *SAconvg,
   int        *NFills,
   double     *Sample
   )
{
  int i, j, rc, nb, Fill;
  double *optPmatp, *optQmatp, *samplep, *dptr;

  double *LEstsptr, *FMatptr;
  int *iptr;
  int *Modelsptr;
  int **BaseModel;

  // main 
  int    MiterP, MiterQ;
  double MxminP, MxminQ, MfminP, MfminQ;

  // --------------------------------------------------------------
  double **sMat = NULL;
  double **iMat = NULL;
  double **dMat = NULL, **indMat = NULL;
  double MaximumValue, MinimumValue;

  LogisticS *logS;

  Pstruct *Pptr;
  Qstruct *Qptr; 

  // --------------------------------------------------------------
  iMat    = mkMatd(*N0, *NT);
  dMat    = mkMatd(*N0, *NT);
  indMat  = mkMatd(*N0, *NT);

  dptr = Mat;
  if ( isnan(*dptr) ) return;

  dptr = Mat;
  MaximumValue = *dptr;
  MinimumValue = *dptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      //  for ( i = 0; i < *N0; i++ ) { 
      //    for ( j = 0; j < *NT; j++ ) {
      dMat[i][j] = *dptr++;
      if ( !isnan(dMat[i][j]) && dMat[i][j] > MaximumValue ) MaximumValue = dMat[i][j];
      if ( !isnan(dMat[i][j]) && dMat[i][j] < MinimumValue ) MinimumValue = dMat[i][j];
      indMat[i][j] = dMat[i][j];
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
                  count, uvals, smallinc);

  rc = load_IF(*N0,*NT,count, uvals, 0, NULL );
  // --------------------------------------------------------------
  
  seed_sgen(*seed1);
  
  // 0.0  Fill in intermittent missing data
  // --------------------------------------------------------------
  
  nb = *nmodel;
  BaseModel = mkMati(*NT,nb);
  
  iptr = inmodel;
  for ( j = 0; j <  nb; j++ ) {
    for ( i = 0; i < *NT; i++ ) { 
      //  for ( i = 0; i < *NT; i++ ) { 
      //    for ( j = 0; j <  nb; j++ ) {
      BaseModel[i][j] = *iptr++;
    }
  }

  logS  = initLogisticS( indMat, *N0, *NT, nb, 25, 10E-10, 10E-10 );

  // for outputs
  LEstsptr  = LEstsM;
  Modelsptr = ModelsM;
  FMatptr   = FMat;
  optPmatp  = outPmatM;
  optQmatp  = outQmatM;

  for ( Fill = 0; Fill < *NFills; Fill++ ) {

    cpMatd( logS->Data,   indMat,    *N0, *NT ); 
    cpMati( logS->Models, BaseModel, *NT,  nb ); 

    rc = FillIn( nb, logS );

    rc = logStoOut( logS, &LEstsptr, &Modelsptr, *NT, nb, 0, Fill+1 );

    cpMatd(iMat,logS->Data,*N0,*NT); 
    for ( i = 0; i < *N0; i++ ) { 
      *FMatptr++ = 0;
      *FMatptr++ = Fill+1;
      for ( j = 0; j < *NT; j++ )
  	*FMatptr++ = iMat[i][j];
    }

    // 1.0
    // --------------------------------------------------------------
    
    // P opt
    rc =  Popt(iMat, &MiterP, &MxminP, &MfminP);
    rc =  toOut( &optPmatp, 0, Fill+1, rc, MiterP, MxminP, MfminP);

    // Q opt
    rc =  Qopt(iMat, &MiterQ, &MxminQ, &MfminQ);
    rc =  toOut( &optQmatp, 0, Fill+1, rc, MiterQ, MxminQ, MfminQ);

    Pptr = IFscr.Pptr;
    rc   = Pinit1(Pptr, iMat, *N0, *NT, 0, 0, 1);
    rc   = updateP( Pptr, MxminP );

    Qptr = IFscr.Qptr;
    rc   = Qinit1(Qptr, iMat, *N0, *NT, 0, 0, 1);
    rc   = updateQ( Qptr, MxminQ );
  }

  // 2.0 -- use last created iMat to generate samples
  // --------------------------------------------------------------

  if ( *NSamples > 0 ) {
    seed_sgen(SEnv.seed0);

    // note NSamles is number of samples here, not number of datasets 
    sMat = mkMatd(*NSamples, *NT);
    rc = Gen_fun( iMat, *N0, *NT, sMat, *NSamples, MxminP , MxminQ );

    // Put the sample into logS and introduce intermittent missing data
    rc = freeMatd(logS->Data);
    rc = freeMati(logS->Last);
    logS->N0   = *NSamples;
    logS->Data = mkMatd( *NSamples, *NT);
    logS->Last = mkMati( *NSamples,   1);

    for ( i = 0; i < *NSamples; i++ ) {
      for ( j = 0; j < *NT; j++ ) {
    	(logS->Data)[i][j] = sMat[i][j];
    	if ( !isnan(sMat[i][j]) ) logS->Last[i][0] = j;
      }
    } 
    mkIM( logS );

    // Write sample
    samplep = Sample;
    for ( i = 0; i < *NSamples; i++ ) {
      *samplep++ = *NFills;
      for ( j = 0; j < *NT; j++ ) {
    	sMat[i][j] = (logS->Data)[i][j];
    	*samplep++ = sMat[i][j];
      }
    }
  }
  rc = distructLogisticS( logS );

  // --------------------------------------------------------------
  freeMatd(iMat    );
  freeMatd(sMat    );
  freeMatd(dMat    );
  freeMatd(indMat  );

  freeMati(BaseModel);

  rc = free_IF();
  rc = free_env();
  freeMatd(uvals);
  return;
}
