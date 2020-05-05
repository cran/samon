// ----------------------------------------------------------------------------
// samon_boot_jk2: Entry point for samon in R
// ----------------------------------------------------------------------------
// 0. Initialize environment
// 1. Popt, Qopt and IF for input data 
// 2. Generate NSamples ( if *NSamples > 0 )
// 3. Popt, Qopt and IF for each bootstrap sample ( if *NSamples > 0 )
// 4. Reconfigure N0 -> N0 - 1 in the environment ( if jk are needed )
// 5. compute jk for input data ( if *Mjk == 1 )
// 6. compute jk for NSamples   ( if *NSamples > 0 and *Sjk == 1 )
// ------------------------------------------------------------------------
// The call:
// All arrays are passed as 1-dim arrays to satisfy R calling convention
// void samon_boot_jk2 ( 
//  int               *N0,   Number of observations
//  int               *NT,   Number of time-points
//  double           *Mat,   *N0 by *NT input matrix
//
//  double      *outPmatM,   Output matrix for main P smoothing parameter
//  double      *outQmatM,   Output matrix for main Q smoothing parameter
//  double     *alphamatM,   Output matrix of IF estimates for main data
//  double    *outPmatMjk,   Repeat for jackknife estimates 
//  double    *outQmatMjk,
//  double   *alphamatMjk,
//  double      *outPmatS,   Output matrix for bootstrap P smoothing
//  double      *outQmatS,   Output matrix for bootstrap Q smoothing
//  double     *alphamatS,   Output matrix of IF estimates for bootstraps
//  double    *outPmatSjk,   Repeat for jackknife estimates of bootstraps
//  double    *outQmatSjk,
//  double   *alphamatSjk,
//
//  int            *seed0,   Seed for generating bootstraps
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
//  int         *NSamples,   Number of bootstrap samples to create
//
//  int          *MaxIter,   Convergence criteria for optimization
//  double       *FAconvg,
//  double       *FRconvg,
//  double       *SAconvg,
//
//  int           *Nalpha,   Number of alphas in alphaList
//  double     *alphaList,   List of alphas
//
//  int              *Mjk,   produce jackknife estimates for main data  
//  int              *Sjk,   produce jackknife estimates for bootstraps
//  int          *retIFiM,   return individual IF values for main
//  int          *retIFiS,   return individual IF valuse for bootstraps
//  double       *IFvalsM,   array to return individual IF values for
//                           main data (if requested)
//  double       *IFvalsS,   array to return individual IF values for
//                           bootstrap data (if requested) 
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
#include "samon_Util.h"
#include "unique.h"

samonEnv SEnv;
IFscratch IFscr;

// ------------------------------------------------------------------------
void samon_boot_jk2 ( 
  int               *N0,   int               *NT,    double          *Mat, 

  double      *outPmatM,   double      *outQmatM,   double     *alphamatM,
  double    *outPmatMjk,   double    *outQmatMjk,   double   *alphamatMjk,
  double      *outPmatS,   double      *outQmatS,   double     *alphamatS,
  double    *outPmatSjk,   double    *outQmatSjk,   double   *alphamatSjk,

  int            *seed0,
  double        *startp,   double         *HSigp,
  double        *startq,   double         *HSigq,
  double            *lb,   double            *ub,
  double         *zeta1,   double         *zeta2,
  int           *NParts,
  int         *NSamples,
  int          *MaxIter,
  double       *FAconvg,   double       *FRconvg,   double       *SAconvg,
  int           *Nalpha,   double     *alphaList,
  int              *Mjk,   int              *Sjk,
  int          *retIFiM,   int          *retIFiS,
  double       *IFvalsM,   double       *IFvalsS,
  int          *lenTfun,   double          *Tfun
                )
{
  int rc, i, j, k ;
  int i0, cnt, bt, btstar;
  double *ptr;

  double *optr;
  double **ifiptr;
  double *optPmatp, *optQmatp;

  double **TfunMat = NULL;

  // main 
  int    MiterP, MiterQ, *MiterPp, *MiterQp;
  double MxminP, MxminQ, *MxminPp, *MxminQp;
  double MfminP, MfminQ, *MfminPp, *MfminQp;

  // bootstrap samples
  int    S1iterP, S1iterQ, *S1iterPp, *S1iterQp;
  double S1xminP, S1xminQ, *S1xminPp, *S1xminQp;
  double S1fminP, S1fminQ, *S1fminPp, *S1fminQp;

  // --------------------------------------------------------------
  double **dMat;
  double **sMat = NULL;
  double **iTmp; 
  double **iSamp;
  double **sample1;

  double MaximumValue, MinimumValue;

  // --------------------------------------------------------------
  iTmp    = mkMatd(*N0, *NT); 
  iSamp   = mkMatd(*N0, *NT);
  sample1 = mkMatd(*N0 - 1, *NT);

  dMat    = mkMatd(*N0, *NT);

  ptr = Mat;
  if ( isnan(*ptr) ) return;

  MaximumValue = *ptr;
  MinimumValue = *ptr;
  for ( j = 0; j < *NT; j++ ) {
    for ( i = 0; i < *N0; i++ ) { 
      //  for ( i = 0; i < *N0; i++ ) { 
      //  for ( j = 0; j < *NT; j++ ) {
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
  optr     = alphamatM;
  ifiptr   = &IFvalsM;
  rc = IF_fun(dMat, MxminP, MxminQ, 0, 0, &optr, *Nalpha, alphaList, *retIFiM, ifiptr);

  // 2.0  The samples
  // --------------------------------------------------------------
  if ( *NSamples > 0 ) { 

    // Generate
    seed_sgen(SEnv.seed0);
 
    sMat = mkMatd(*NSamples * *N0, *NT);
    rc = Gen_fun( dMat, *N0, *NT, sMat, *NSamples * *N0, MxminP , MxminQ );

  // --------------------------------------------------------------

    S1iterPp = &S1iterP;
    S1xminPp = &S1xminP;
    S1fminPp = &S1fminP; 
    S1iterQp = &S1iterQ;
    S1xminQp = &S1xminQ;
    S1fminQp = &S1fminQ; 

    optr     = alphamatS;
    optPmatp = outPmatS;
    optQmatp = outQmatS;

    // 3.0 IF for each bootstrap sample
    //     --------------------------------
    for ( bt = 0; bt < *NSamples; bt++ ) {

      // copy sample for IF computation
      btstar = bt * *N0; 
      i0 = 0;
      for ( i = btstar; i < btstar + *N0; i++ ) {
	for ( j = 0; j < *NT; j++ ) 
  	  iSamp[i0][j] = sMat[i][j];
	i0++;
      }

      // P opt
      rc =  Popt(iSamp, S1iterPp, S1xminPp, S1fminPp);
      rc =  toOut( &optPmatp, bt+1, 0, rc, S1iterP, S1xminP, S1fminP);
      SEnv.minSampPQ[bt][0] = S1xminP;

      // Q optimization
      rc =  Qopt(iSamp, S1iterQp, S1xminQp, S1fminQp);
      rc =  toOut( &optQmatp, bt+1, 0, rc, S1iterQ, S1xminQ, S1fminQ);
      SEnv.minSampPQ[bt][1] = S1xminQ;

      // IF estimates
      ifiptr   = &IFvalsS;
      rc = IF_fun(iSamp, S1xminP, S1xminQ, 0, bt+1, 
  		  &optr, *Nalpha, alphaList, *retIFiS, ifiptr);
    }
  }

  // 4.0 Restructure
  // --------------------------------------------------------------
  if ( *Mjk == 1 || *Sjk == 1 ) {
    rc = reload_env ( *N0 - 1, MxminP, MxminQ );
    rc = free_IF();
    rc = load_IF(*N0 - 1,*NT, count, uvals, *lenTfun, TfunMat) ;
  }
  // --------------------------------------------------------------

  // 5.0 Jackknife original (main) data
  if ( *Mjk == 1 ) {

    optr     = alphamatMjk;
    optPmatp = outPmatMjk;
    optQmatp = outQmatMjk;

    if ( SEnv.SmallV < SEnv.minMainPQ[0] && SEnv.minMainPQ[0] < SEnv.HSigp ) SEnv.startp = SEnv.minMainPQ[0];
    if ( SEnv.SmallV < SEnv.minMainPQ[1] && SEnv.minMainPQ[1] < SEnv.HSigq ) SEnv.startq = SEnv.minMainPQ[1];

    for ( i = 0; i < *N0; i++ ) {
      cnt = 0;
      for ( j = 0; j < *N0; j++ ) {
  	if ( j != i ) {
  	  for ( k = 0; k < *NT; k++ ) 
  	    sample1[cnt][k] = dMat[j][k];
  	  cnt++;
  	}
      } 

      // optimize P
      rc =  Popt(sample1, S1iterPp, S1xminPp, S1fminPp);
      rc =  toOut( &optPmatp, 0, i+1, rc, S1iterP, S1xminP, S1fminP);

      // optimize Q
      rc =  Qopt(sample1, S1iterQp, S1xminQp, S1fminQp);
      rc =  toOut( &optQmatp, 0, i+1, rc, S1iterQ, S1xminQ, S1fminQ);

      // IF estimates
      rc = IF_fun( sample1, S1xminP, S1xminQ, i+1, 0, &optr, *Nalpha, alphaList, 0, 0);
    }
  }
  // --------------------------------------------------------------
  if ( *Sjk == 1 ) {

    optr     = alphamatSjk;
    optPmatp = outPmatSjk;
    optQmatp = outQmatSjk;

    // 6.0 Jackknife for each bootstrap sample
    // ---------------------------------------
    for ( bt = 0; bt < *NSamples; bt++ ) {

      if ( SEnv.SmallV < SEnv.minSampPQ[bt][0] && SEnv.minSampPQ[bt][0] < SEnv.HSigp )
  	SEnv.startp = SEnv.minSampPQ[bt][0];
      else if ( SEnv.SmallV < SEnv.minMainPQ[0] && SEnv.minMainPQ[0] < SEnv.HSigp )
  	SEnv.startp = SEnv.minMainPQ[0];
      else SEnv.startp = *startp;

      if ( SEnv.SmallV < SEnv.minSampPQ[bt][1] && SEnv.minSampPQ[bt][1] < SEnv.HSigq )
  	SEnv.startq = SEnv.minSampPQ[bt][1];
      else  if ( SEnv.SmallV < SEnv.minMainPQ[1] && SEnv.minMainPQ[1] < SEnv.HSigq )
  	SEnv.startq = SEnv.minMainPQ[1];
      else SEnv.startq = SEnv.minMainPQ[1];

      // copy sample for jackknife
      btstar = bt * *N0; 
      i0 = 0;
      for ( i = btstar; i < btstar + *N0; i++ ) {
  	for ( j = 0; j < *NT; j++ ) 
  	  iTmp[i0][j] = sMat[i][j];
  	i0++;
      }

      for ( i = 0; i < *N0; i++ ) {

  	cnt = 0;
  	for ( j = 0; j < *N0; j++ ) {
  	  if ( j != i ) {
  	    for ( k = 0; k < *NT; k++ ) {
  	      sample1[cnt][k] = iTmp[j][k];
  	    }
  	    cnt++;
  	  }
  	} 

	// P opt
  	rc =  Popt(sample1, S1iterPp, S1xminPp, S1fminPp);
  	rc =  toOut( &optPmatp, bt+1, i+1, rc, S1iterP, S1xminP, S1fminP);

	// Q opt
  	rc =  Qopt(sample1, S1iterQp, S1xminQp, S1fminQp);
  	rc =  toOut( &optQmatp, bt+1, i+1, rc, S1iterQ, S1xminQ, S1fminQ);

	// IF
  	rc = IF_fun( sample1, S1xminP, S1xminQ, i+1, bt+1, &optr, *Nalpha, alphaList, 0, 0);
      }
    }
  }

  // --------------------------------------------------------------

  freeMatd( dMat );
  freeMatd( iTmp    );
  freeMatd( iSamp   );
  freeMatd( sample1 );

  if ( *NSamples > 0 ) {
    freeMatd(sMat);
  }

  rc = free_IF();
  rc = free_env();
  freeMatd( uvals );
  if ( *lenTfun > 0 ) freeMatd( TfunMat );
  return;
}
// --------------------------------------------------------------
