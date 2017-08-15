// ----------------------------------------------------------------------------
// samon_ngenIMIF: Entry point for samon in R -- intermittent missing data
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
//  double         *FMatM,   Filled in main data
//  double        *LEstsM,   Logistic results for main data
//  int          *ModelsM,   Models for main data
//
//  double         *FMatS,   Filled in sample data
//  double        *LEstsS,   Logistic results for sample data
//  int          *ModelsS,   Models for sample data
//
//  double      *outPmatM,   Output matrix for main P smoothing parameter
//  double      *outQmatM,   Output matrix for main Q smoothing parameter
//  double     *alphamatM,   Output matrix of IF estimates for main data
//
//  double      *outPmatS,   Output matrix for sample P smoothing parameter
//  double      *outQmatS,   Output matrix for sample Q smoothing parameter
//  double     *alphamatS,   Output matrix of IF estimates for sample data
//
//  int            *seed0,   Seed for generating bootstraps (not used here)
//  int            *seed1,   Seed filling in intermittent missing data
//  double        *startp,   Start value to use for P smoothing parameter
//  double         *HSigp,   High value for P smoothing
//  double        *startq,   Start value tor use for Q smoothing parameter
//  double         *HSigq,   High value for Q smoothing
//
//  double            *lb,   parameters for cumulative beta function
//  double            *ub,
//  double         *zeta1,
//  double         *zeta2,
//
//  int           *NParts,   Number of partitions
//  int         *NSamples,   Number of samples to make
//
//  int          *MaxIter,   Convergence criteria for optimization
//  double       *FAconvg,
//  double       *FRconvg,
//  double       *SAconvg,
//
//  int           *Nalpha,   Number of alphas in alphaList
//  double     *alphaList,   List of alphas
//
//  int           *NFills,   Number of fillins to do
//
//  int          *retIFiM,   logical return individual IF values
//  int          *retIFiS,   logical return individual IF values (samples)
//  double       *IFvalsM,   individual IF values for (if requested)
//  double       *IFvalsS,   individual IF values for (if requested)
//  
//  int        *retSample,   Should the sample be returned
//  double        *Sample,   The output sample
//
//  int         *retFMatM,   Should FMatM be returned
//  int         *retFMatS,   Should FMatS be returned
//
//  int          *lenTfun,   Tilt function nrow
//  double          *Tfun,   Tilt matrix (*lenTfun by 2)
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
void samon_ngenIMIF ( 
   int            *N0,  int           *NT,  double       *Mat, 
   int        *nmodel,  int      *inmodel,
   double      *FMatM,  double    *LEstsM,  int      *ModelsM,
   double      *FMatS,  double    *LEstsS,  int      *ModelsS,
   double   *outPmatM,  double  *outQmatM,  double *alphamatM,
   double   *outPmatS,  double  *outQmatS,  double *alphamatS,
   int         *seed0,  int        *seed1,
   double     *startp,  double     *HSigp,
   double     *startq,  double     *HSigq,
   double         *lb,  double        *ub,
   double      *zeta1,  double     *zeta2,
   int        *NParts,
   int      *NSamples, 
   int       *MaxIter,  double   *FAconvg,  double   *FRconvg,  double  *SAconvg,
   int        *Nalpha,  double *alphaList,
   int        *NFills,
   int       *retIFiM,  int      *retIFiS,
   double    *IFvalsM,  double   *IFvalsS,
   int     *retSample,
   double     *Sample,
   int      *retFMatM,  int     *retFMatS, 
   int       *lenTfun,  double      *Tfun
   )
{
  int i, j, rc, nb, Fill, samp;
  double *ptr, *optr, *optPmatp, *optQmatp, *samplep, *dptr;

  double **ifiptr;

  double *LEstsptr, *FMatptr;
  int *iptr;
  int *Modelsptr;
  int **BaseModel;

  double **TfunMat = NULL;

  // main 
  int    MiterP, MiterQ;
  double MxminP, MxminQ, MfminP, MfminQ;

  // --------------------------------------------------------------
  double **sMat = NULL;
  double **iMat = NULL;
  double **dMat = NULL, **indMat = NULL;
  double MaximumValue, MinimumValue;

  LogisticS *logS;

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
    
  if ( *lenTfun > 0 ) {
    TfunMat = mkMatd(*lenTfun,2);
    ptr = Tfun;
    for ( j = 0; j < 2; j++ ) {
      for ( i = 0; i < *lenTfun; i++ ) {
    //    for ( i = 0; i < *lenTfun; i++ ) {
    //      for ( j = 0; j < 2; j++ ) {
	TfunMat[i][j] = *ptr++;
      }
    }
  }

  rc = load_env ( *N0, *NT, *seed0, *startp, *HSigp, *startq, *HSigq,
                  *lb, *ub, *zeta1, *zeta2,
                  *NParts, *NSamples,
                  *MaxIter, *FAconvg, *FRconvg, *SAconvg, 
                  count, uvals, smallinc/4.0);

  rc = load_IF(*N0,*NT,count, uvals, *lenTfun, TfunMat );
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
  FMatptr   = FMatM;
  optPmatp  = outPmatM;
  optQmatp  = outQmatM;
  optr      = alphamatM;
  ifiptr    = &IFvalsM;


  for ( Fill = 0; Fill < *NFills; Fill++ ) {

    cpMatd( logS->Data,   indMat,    *N0, *NT ); 
    cpMati( logS->Models, BaseModel, *NT,  nb ); 

    rc = FillIn( nb, logS );

    rc = logStoOut( logS, &LEstsptr, &Modelsptr, *NT, nb, 0, Fill+1 );

    cpMatd(iMat,logS->Data,*N0,*NT); 
    if ( *retFMatM == 1 ) { 
      for ( i = 0; i < *N0; i++ ) { 
	*FMatptr++ = 0;
	*FMatptr++ = Fill+1;
	for ( j = 0; j < *NT; j++ )
	  *FMatptr++ = iMat[i][j];
      }
    }

    // 1.0
    // --------------------------------------------------------------
    
    // P opt
    rc =  Popt(iMat, &MiterP, &MxminP, &MfminP);
    rc =  toOut( &optPmatp, 0, Fill+1, rc, MiterP, MxminP, MfminP);

    // Q opt
    rc =  Qopt(iMat, &MiterQ, &MxminQ, &MfminQ);
    rc =  toOut( &optQmatp, 0, Fill+1, rc, MiterQ, MxminQ, MfminQ);

    // IF
    rc = IF_fun(iMat, MxminP, MxminQ, Fill+1, 0, &optr, *Nalpha, alphaList, *retIFiM, ifiptr);
  }

  // 2.0 -- use last created iMat to generate samples
  // --------------------------------------------------------------

  if ( *NSamples > 0 ) {
    seed_sgen(SEnv.seed0);

    sMat = mkMatd(*NSamples * *N0, *NT);
    rc = Gen_fun( iMat, *N0, *NT, sMat, *NSamples * *N0, MxminP , MxminQ );

    // Put the sample into logS and introduce intermittent missing data
    rc = freeMatd(logS->Data);
    rc = freeMati(logS->Last);
    logS->N0   = *NSamples * *N0;
    logS->Data = mkMatd( *NSamples * *N0, *NT);
    logS->Last = mkMati( *NSamples * *N0,   1);

    for ( i = 0; i < *NSamples * *N0; i++ ) {
      for ( j = 0; j < *NT; j++ ) {
    	(logS->Data)[i][j] = sMat[i][j];
    	if ( !isnan(sMat[i][j]) ) logS->Last[i][0] = j;
      }
    } 
    mkIM( logS );

    // Write sample
    samplep = Sample;
    for ( i = 0; i < *NSamples * *N0; i++ ) {
      if ( *retSample == 1 ) *samplep++ = floor(i / *N0) + 1;
      for ( j = 0; j < *NT; j++ ) {
	sMat[i][j] = (logS->Data)[i][j];
	if ( *retSample == 1 ) *samplep++ = sMat[i][j];
      }
    }
  }
  rc = distructLogisticS( logS );

  if ( *NSamples > 0 ) {

    // 3.0 For each sample fill n times and Popt, Qopt and IF
    // --------------------------------------------------------------
    // outputs
    LEstsptr  = LEstsS;
    Modelsptr = ModelsS;
    FMatptr   = FMatS;
    optPmatp  = outPmatS;
    optQmatp  = outQmatS;
    optr      = alphamatS;
    ifiptr    = &IFvalsS;

    logS  = initLogisticS( indMat, *N0, *NT, nb, 25, 10E-10, 10E-10 );

    for ( samp=0; samp < *NSamples; samp++ ) {

      // read sample for processing
      for ( i = 0; i < *N0; i++ ) {
        for ( j = 0; j < *NT; j++ ) {
	  dMat[i][j] = sMat[i+(samp * *N0)][j];
	  if ( !isnan(dMat[i][j]) ) logS->Last[i][0] = j;
        }
      } 

      for ( Fill = 0; Fill < *NFills; Fill++ ) {

        cpMatd( logS->Data,   dMat,      *N0, *NT ); 
        cpMati( logS->Models, BaseModel, *NT,  nb );  // for now

        rc = FillIn( nb, logS );

        rc = logStoOut( logS, &LEstsptr, &Modelsptr, *NT, nb, samp+1, Fill+1 );

        cpMatd(iMat,logS->Data,*N0,*NT); 
	if ( *retFMatS == 1 ) {
	  for ( i = 0; i < *N0; i++ ) { 
	    *FMatptr++ = samp+1;
	    *FMatptr++ = Fill+1;
	    for ( j = 0; j < *NT; j++ )
	      *FMatptr++ = iMat[i][j];
	  }
	}

        // 3.0 b
        // --------------------------------------------------------------
        // P opt
        rc =  Popt(iMat, &MiterP, &MxminP, &MfminP);
        rc =  toOut( &optPmatp, samp+1, Fill+1, rc, MiterP, MxminP, MfminP);

        // Q opt
        rc =  Qopt(iMat, &MiterQ, &MxminQ, &MfminQ);
        rc =  toOut( &optQmatp, samp+1, Fill+1, rc, MiterQ, MxminQ, MfminQ);

        // IF
        rc = IF_fun(iMat, MxminP, MxminQ, Fill+1, samp+1, &optr, *Nalpha, alphaList, *retIFiS, ifiptr);
      }
    }
    rc = distructLogisticS( logS );
  }
  // --------------------------------------------------------------
  freeMatd(iMat    );
  freeMatd(sMat    );
  freeMatd(dMat    );
  freeMatd(indMat  );

  freeMati(BaseModel);

  rc = free_IF();
  rc = free_env();
  freeMatd(uvals);

  if ( *lenTfun > 0 ) freeMatd( TfunMat );
  return;
}
