// ----------------------------------------------------------------------------
// IF: computes influence function estimates.
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"
#include "Qstruct.h"
#include "IFscratch.h"
#include "IF_fun_beta.h"
#include "beta_cdf.h"

double **tPMat;
double **tdv;

double **tInt;
double **tBInt;
double **tCInt; 

double **tPre;

double **tT2;
double **tB3_3;
double **tCM;

int MXV;
int MXT;

int N0, NT;

// ------------------------------------------------------------------------

double UU;

double *tEalpha;
double *trfun;

int **tV;
double **tYY0;

double *tQ0;
double **Top;
double **Bot;

int position( double *arr, int n, double Y );

int IF_fun(double **y, double in_sigmap, double in_sigmaq, int tp1, int rep, double **OOptr, int Nalpha, double *alphalist, int ifi,  double **ifiptr)
{
  int i, k;
  int Firstalpha;
  double alpha; 

  double UU, VarUU;
  double  A, VarA;

  double Sigmap, Sigmaq;

  N0 = SEnv.N0;
  NT = SEnv.NT;

// --------------------------------------------------------------------

  Pstruct *Pptr;
  Qstruct *Qptr; 

  Pptr = IFscr.Pptr;
  Pinit1(Pptr, y, N0, NT, 0, 0, 1);
  updateP( Pptr, in_sigmap );

  Qptr = IFscr.Qptr;
  Qinit1(Qptr, y, N0, NT, 0, 0, 1);
  updateQ( Qptr, in_sigmaq );

// --------------------------------------------------------------------

  MXV = SEnv.nunique;
  MXT = NT+1;

   tPMat   = IFscr.tPMat;
   tdv     = IFscr.tdv;
   tInt    = IFscr.tInt;
   tBInt   = IFscr.tBInt;
   tCInt   = IFscr.tCInt;
   tPre    = IFscr.tPre;
   tT2     = IFscr.tT2;
   tB3_3   = IFscr.tB3_3;
   tCM     = IFscr.tCM;
   Top     = IFscr.Top;
   Bot     = IFscr.Bot;

   tV      = IFscr.tV;
   zeroi( tV,    MXV, MXT );

   tYY0    = y;   // copy pointer 

   tQ0     = IFscr.tQ0;
   tEalpha = IFscr.tEalpha;
   trfun   = IFscr.trfun; 

   for ( i = 0; i < MXV; i++ ) {
     tQ0[i]     = 0.0;
     tEalpha[i] = 0.0;
   }

   // Q0
   // -----------------------------------------------------------------
   mkQ0(Pptr);

   Sigmap = in_sigmap;
   Sigmaq = in_sigmaq;

   Firstalpha = 1;
   for ( k = 0; k < Nalpha; k++ ) {

     alpha = alphalist[k];

     // Ealpha 
     // ------------------------
     for ( i=0; i < MXV; i++ ) {
       tEalpha[i] = exp( alpha * trfun[i]);
     }

     if ( Firstalpha == 1 ) {
       mkPQMat( Pptr, Qptr, MXV, MXV, NT, alpha );
       updateT( Qptr, alpha );
       Firstalpha = 0;
     } else {
       updateT( Qptr, alpha );
     }

     mkInt(Pptr,Qptr);
     mkPre(Pptr,Qptr);
     mkH(Qptr); 
     mkT2();
     mkB3();
     mkBMat(Qptr);
     mkC();

     mkUU(rep, Sigmap, Sigmaq, alpha, &UU, &VarUU, &A, &VarA, Qptr, ifi, ifiptr);

     // Write to return matrix  
     *((*OOptr)++)  = rep;
     *((*OOptr)++)  = tp1;
     *((*OOptr)++)  = alpha;
     *((*OOptr)++)  = A;
     *((*OOptr)++)  = VarA;
     *((*OOptr)++)  = UU;
     *((*OOptr)++)  = VarUU;
   }
   return 0;
} 

// --------------------------------------------------------------------
 int mkQ0(Pstruct *X)
 {
   int i, nb, ip;

   // ----- Calculate Q0;
   for ( i=0; i < MXV; i++ ) {
     tQ0[i] = 0.0;
   }

   nb = (X->Nb)[0];
   for ( i=0; i < nb; i++ ) {
     ip = ((X->Posb)[0])[i];   // position  - integer
     tQ0[ ip ] = ( ((X->b)[0])[i][1] + ((X->b)[0])[i][2] ) / N0;
   }
   return 0;
 }

 // --------------------------------------------------------------------

 int mkPQMat( Pstruct *Pptr, Qstruct *Qptr, int nr, int nc, int nt, double alpha )
 {
   int i,t, ip;

   for ( i = 0; i < nr; i++ ) { 
     for ( t = 0; t < (nt); t++ ) {
       tPMat[i][ t] = 0.0;
     }
   }

   for ( t = 0; t < (nt-1); t++ ) {
     for ( i = 0; i < (Pptr->Nb)[t]; i++ ) {
       ip = ((Pptr->Posb)[t])[i];
       tPMat[ ip ][t] = ((Pptr->P)[t])[i];
     }
   }
   return 0;
 }

 // --------------------------------------------------------------------

 int updateT(Qstruct* Qptr, double alpha)
 {
   int t, i, j;
   double sst;
   int nr, nc, nt;
   int ip, jp;
   double **Q, **TQ;

   nt = Qptr->NT;

   for ( t=0; t < (nt-1); t++ ) {
     for ( i = 0; i < MXV; i ++ ) {
       tdv[i][t] = 0.0;
     }
   }

   // tilted by alpha
   for ( t=0; t < (nt-1); t++ ) {
     nr = (Qptr->Nr)[t];
     nc = (Qptr->Nc)[t];
     Q  = (Qptr->Q )[t];
     TQ = (Qptr->TQ)[t];
     for ( i=0; i < nr; i++ ) {
       ip = ((Qptr->Posr)[t])[i];
       sst = 0.0;
       for ( j=0; j < nc; j++ ) {
	 jp = ((Qptr->Posc)[t])[j];
	 TQ[i][j] = Q[i][j] * exp ( alpha * trfun[ jp ] );
	 sst = sst + TQ[i][j];
       }
       if ( sst > 0.0 ) {
	 tdv[ip][t] = sst;
	 for ( j=0; j < nc; j++ ) {
	   TQ[i][j] = TQ[i][j] / sst;
	 }
       }
     }
   }  
   return 0;
 }

 // --------------------------------------------------------------------

 int mkInt( Pstruct *Pptr, Qstruct *Qptr)
 {
   int i, j, t;
   double iv;
   int ip, jp;
   int nr, nc;
   double **Q, **TQ;

   // --------------------------------------------------------------------
   for ( t = 0; t < NT; t++ ) {
     for (i = 0; i < MXV; i++ ) {
       tInt[i][t]   = 0;
       tBInt[i][t]  = 0;
       tCInt[i][t]  = 0;
     }
   }

   for( t = (NT-1); t > -1; t-- ) {

     if ( t == (NT-1) ) {
       for ( i = 0; i < (Qptr->Nc)[t-1]; i++ ) {
	 iv = (Qptr->Qc)[t-1][i];  
	 ip = (Qptr->Posc)[t-1][i];  
	 tInt[ ip ][t] = iv;
       }
     } else {
       nr = (Qptr->Nr)[t];
       nc = (Qptr->Nc)[t];
       Q  = (Qptr->Q )[t];
       TQ = (Qptr->TQ)[t]; 

       for ( i = 0; i < nr; i++ ) {
	 iv  = (Qptr->Qr)[t][i];  
	 ip  = (Qptr->Posr)[t][i];  

	   for ( j = 0; j < nc; j++ ) {
	     jp    = (Qptr->Posc)[t][j];

	     if ( Q[i][j] > 0.0 ) {

	       tInt[ip][t] = tInt[ip][t] +
		 tPMat[ip][t]   * Q [i][j] * tInt[jp][t+1] +
		 ( 1 - tPMat[ip][t] ) * TQ[i][j] * tInt[jp][t+1];

	       if ( tdv[ip][t] != 0 )
		 tBInt[ip][t] = tBInt[ip][t] +
		   ( 1 - tPMat[ip][t] ) * Q[i][j] * ( 1 / tdv[ip][t] ) * ( 1/ tdv[ip][t] ) * tEalpha[jp] * tInt[jp][t+1];

	       if ( tdv[ip][t] != 0 )
		 tCInt[ip][t] = tCInt[ip][t] +
		   Q[i][j] * tInt[jp][t+1] -
		   Q[i][j] * ( 1 / tdv[ip][t] ) * tEalpha[jp] * tInt[jp][t+1];
	     }
	   }
       }
     }
   }
   return 0;
 } 

 // --------------------------------------------------------------------

 int mkPre(Pstruct *Pptr, Qstruct *Qptr) 
 {
   int i, j, t;
   int nr, nc;
   int ip, jp;
   double **Q, **TQ;

  // --------------------------------------------------------------------

   for( t = 0; t <= NT; t++ ) {
     for ( i = 0; i < MXV; i++ ) {
       if ( t == 0 ) {
	 Top[i][t] = tQ0[i]; 
	 Bot[i][t] = tQ0[i];
	 tPre[i][t] = 1;
       }
       else {
	 Top[i][t] = 0;
	 Bot[i][t] = 0;
	 tPre[i][t] = 0;
       }
     }
   }

   for( t = 1; t <= NT; t++ ) {
     nr  = (Qptr->Nr)[t-1];
     nc  = (Qptr->Nc)[t-1];
     Q   = (Qptr->Q )[t-1];
     TQ  = (Qptr->TQ)[t-1]; 

     for (i = 0; i < nc; i++ ) {
       ip  = (Qptr->Posc)[t-1][i];  

	 for (j = 0; j < nr; j++ ) {

	   jp  = (Qptr->Posr)[t-1][j];  
	   if ( Q[j][i] > 0.0 ) {
	     Top[ip][t] = Top[ip][t] +
	       tPMat[jp][t-1] * Q[j][i] * Top[jp][t-1] +
	       ( 1 - tPMat[jp][t-1] ) * TQ[j][i] * Top[jp][t-1];
	     Bot[ip][t] = Bot[ip][t] +
	       tPMat[jp][t-1] * Q[j][i] * Bot[jp][t-1];
	   }
	 }
	 //      }
       if ( Bot[ip][t] != 0 ) tPre[ip][t] =  Top[ip][t] / Bot[ip][t];
     }
   }
   return 0;
 } 

 // --------------------------------------------------------------------

 int mkH(Qstruct *Qptr)
 {
   int i, j, t;
   int nr, nc;
   int ip, jp;
   double **H; 


   for ( t=0; t < (NT-1); t++ ) {
     nr  = (Qptr->Nr)[t];
     nc  = (Qptr->Nc)[t];
     H   = (Qptr->H )[t];

     for ( i = 0; i < nr; i++ ) {
       ip = (Qptr->Posr)[t][i];  

       if ( tPMat[ip][t] != 0 && tdv[ip][t] != 0 ) {
	 for ( j = 0; j < nc; j++ ) {
	   jp = (Qptr->Posc)[t][j];  
	   H[i][j] =  tPre[ip][t] * ( 1 + (( 1 - tPMat[ip][t] ) / (tPMat[ip][t])) * (1/tdv[ip][t]) * tEalpha[jp] ) * ( tInt[jp][t+1] );
	 }
       }
     }
   }
   return 0;
 }

 // --------------------------------------------------------------------
 int mkT2(void)
 {
   int i, t;

   for ( t = 0; t < (NT-1); t++ ) {
     if ( t == 0 ) {
       for ( i = 0; i < MXV; i++ ) {
	 if ( tPMat[i][t] != 0 ) tT2[i][t+1] = tInt[i][t] / tPMat[i][t];
	 else tT2[i][t+1] = 0;
       }
     }
     else {
       for ( i = 0; i < MXV; i++ ) {
	 tT2[i][t+1] = 0;
	 if ( tPMat[i][t] != 0) {
	   tT2[i][t+1] = tPre[i][t] * tInt[i][t] / tPMat[i][t];
	 }
       }
     }
   }
   return 0;
 }
 // --------------------------------------------------------------------

 int mkB3(void)
 {
   int i, t;

   for ( t = 0; t < (NT-1); t++ ) {
     if ( t == 0 ) {
       for ( i = 0; i < MXV; i++ ) {
	 if ( tPMat[i][t] != 0 ) tB3_3[i][t+1] = tBInt[i][t] / tPMat[i][t];
	 else tB3_3[i][t+1] = 0;
       }
     }
     else {
       for ( i = 0; i < MXV; i++ ) {
	 tB3_3[i][t+1] = 0;
	 if ( tPMat[i][t] != 0) {
	   tB3_3[i][t+1] = tPre[i][t] * tBInt[i][t] / tPMat[i][t];
	 }
       }
     }
   }
   return 0;
 }

 // --------------------------------------------------------------------

 int mkC(void)
 {
   int i, t;

   for ( t = 0; t < NT; t++ ) {
     if ( t == 0 ) {
       for ( i = 0; i < MXV; i++ ) {
	 tCM[i][t+1] = tCInt[i][t];
       }
     }
     else {
       for ( i = 0; i < MXV; i++ ) {
	 tCM[i][t+1] = 0;
	 tCM[i][t+1] = tPre[i][t] * tCInt[i][t];
       }
     }
   }
   return 0;
  }

 // --------------------------------------------------------------------

 int mkBMat(Qstruct *Qptr)
 {
   int i, j, t;
   int nr, nc;
   int ip, jp;
   double **H, **IFB;

   for( t=0; t < NT; t++ ) {
     nr  = (Qptr->Nr  )[t];
     nc  = (Qptr->Nc  )[t];
     H   = (Qptr->H   )[t];
     IFB = (Qptr->IFB )[t];

     for (i=0; i < nr; i++ ) {
       ip = (Qptr->Posr)[t][i];  
       for (j=0; j< nc; j++ ) {
	 jp = (Qptr->Posc)[t][j];  
	 IFB[i][j] = H[i][j] - tT2[ip][t+1] - tB3_3[ip][t+1] * tEalpha[jp] + tB3_3[ip][t+1] * tdv[ip][t];
       }
     }
   }
   return 0;
 }

 // --------------------------------------------------------------------

 int mkUU(int rep, double sigmap, double sigmaq, double alpha, double *UUptr, double *UU2ptr, double *Aptr, double *A2ptr, Qstruct *Qptr, int ifi, double **ifiptr )
 {
   int i, t, ip, ipr, ipc;
   double UU;
   double sumUU, sumUU2, meanUU, VarUU;
   double sumA,  sumA2,  meanA,  VarA;

   sumA   = 0.0;
   sumUU  = 0.0;
   sumA2  = 0.0;
   sumUU2 = 0.0;
   for ( i = 0; i < N0; i++ ) {
     // look at tYY0[i][0] and find its location
     ip = position( SEnv.uvalues, SEnv.nunique, tYY0[i][0] );
     UU = tInt[ip][0];
     sumA   = sumA  + UU;
     sumA2  = sumA2 + UU * UU;

     if ( ifi == 1 ) {
       *((*ifiptr)++) =   rep;
       *((*ifiptr)++) = alpha;
       *((*ifiptr)++) =     i;
       *((*ifiptr)++) =    UU;
     }

     for ( t = 0; t < (NT-1); t++ ) {
       if ( !isnan(tYY0[i][t+1]) )  {
	 ipr = position( (Qptr->Qr)[t], (Qptr->Nr)[t], tYY0[i][t] );
	 ipc = position( (Qptr->Qc)[t], (Qptr->Nc)[t], tYY0[i][t+1] );
	 UU = UU + ((Qptr->IFB)[t])[ ipr ][ ipc ];
       }

       if ( !isnan(tYY0[i][t]) ) {
	   ip  = position( SEnv.uvalues, SEnv.nunique, tYY0[i][t] );
	   if ( !isnan(tYY0[i][t+1]) ) UU = UU + ( 1 - tPMat[ ip ][t] ) * tCM[ ip ][t+1] ;   
	   if (  isnan(tYY0[i][t+1]) ) UU = UU + ( 0 - tPMat[ ip ][t] ) * tCM[ ip ][t+1] ;   
       }
     }

     if ( ifi == 1 ) *((*ifiptr)++) =  UU;

     sumUU  = sumUU + UU;
     sumUU2 = sumUU2 + UU * UU;
   }

   meanA  = sumA / N0;
   meanUU = sumUU / N0;
   VarA   = ((sumA2  / N0) - meanA  * meanA  ) / N0;
   VarUU  = ((sumUU2 / N0) - meanUU * meanUU ) / N0;

   *UUptr  = meanUU;
   *UU2ptr = VarUU;
   *Aptr   = meanA;
   *A2ptr  = VarA;
 
  return 0;
}
// ------------------------------------------------------------------------

int position( double *arr, int n, double Y )
{
  int i;

  if ( Y < arr[0]   ) return -1;
  if ( Y > arr[n-1] ) return -1; 
  for ( i = 0; i < n; i++ ) {
    if ( Y == arr[i] ) return i;
  }
  return -1;
}
