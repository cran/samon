// ----------------------------------------------------------------------------
// Create, populate, maintain and delete Pstruct objects
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Pstruct.h"

double lossP( double sigma, int NParts, double *deriv1, double *deriv2 );

// ------------------------------------------------------------------------

// Create a Pstruct but leave the population til later when the data becomes available
Pstruct *Pinit0( int N0, int NT, int  size, int type ) {

  int t,i, MV;
  Pstruct *Xptr;

  Xptr = malloc(sizeof(Pstruct));
  Xptr->NT   = NT;
  Xptr->N    = N0;
  Xptr->Type = type;

  if ( 0 < size && size < N0 ) MV = size;
  else MV = N0;

  Xptr->Na  = malloc( NT * sizeof(int));
  Xptr->Nb  = malloc( NT * sizeof(int));

  for ( i=0; i < NT; i++ ) {
    Xptr->Na[i]  = 0;
    Xptr->Nb[i]  = 0;
  }

  (Xptr->a    ) = malloc( NT * sizeof( double** ));
  (Xptr->b    ) = malloc( NT * sizeof( double** ));
  (Xptr->P    ) = malloc( NT * sizeof( double*  ));
  (Xptr->Posb ) = malloc( NT * sizeof( int*     ));

  if ( type == 0 ) {
    (Xptr->D1) = malloc( NT * sizeof( double* ));
    (Xptr->D2) = malloc( NT * sizeof( double* ));
  }

  for ( t = 0; t < NT-1; t++ ) {

    Xptr->a[t] = mkMatd( MV, 3 );
    Xptr->b[t] = mkMatd( MV, 3 );

    (Xptr->P )[t] = malloc( MV * sizeof( double ));
    if ( type == 0 ) {
      (Xptr->D1)[t] = malloc( MV * sizeof( double ));
      (Xptr->D2)[t] = malloc( MV * sizeof( double ));
    }
    (Xptr->Posb )[t] = malloc( MV * sizeof( int   ));
  }
  return Xptr;
}

// destruct a Pstruct
int Pdestruct( Pstruct *Xptr ) {

  int t, type;
  int NT;

  type = Xptr->Type; 
  NT   = Xptr->NT;

  if ( NT == 0 ) return 0;

  for ( t = 0; t < (NT-1); t++ ) {
    freeMatd( (Xptr->a)[t] );
    freeMatd( (Xptr->b)[t] );

    free( (Xptr->P    )[t] );
    free( (Xptr->Posb )[t] );

    if ( type == 0 ) {
     free( (Xptr->D1)[t] );
     free( (Xptr->D2)[t] );
    }
  }

  free( (Xptr->a    ) ); 
  free( (Xptr->b    ) ); 
  free( (Xptr->P    ) ); 
  free( (Xptr->Posb ) );

  if ( type == 0 ) {
    free( (Xptr->D1) ); 
    free( (Xptr->D2) ); 
  }
  free( (Xptr->Na ) );
  free( (Xptr->Nb ) );

  free(  Xptr       );
  return 0;
}


// populate an existing Pstruct -- i.e. refill the a and b tables 
// the change here is to support non integer Y values
int Pinit1( Pstruct *Xptr, double **Y, int N0, int NT, int start_cut, int stop_cut, int type ) {

  int t,i, k;
  int tcnta, tcntb;
  int na, nb;
  double **zmata, **zmatb;
  int *nptra, *nptrb;

  for ( i=0; i < NT; i++ ) {
    Xptr->Na[i]  = 0;
    Xptr->Nb[i]  = 0;
  }

  // these are now double**
  zmata = (SEnv.Pscrsch)->zmata;
  zmatb = (SEnv.Pscrsch)->zmatb;

  for ( t = 0; t < NT-1; t++ ) {

    tcnta = 0;
    tcntb = 0;
    for ( i = 0; i < N0; i++ ) {
      if ( !isnan(Y[i][t]) ) {
	if ( type == 0 ) {
	  if ( start_cut <= i && i <= stop_cut ) {
	    zmatb[tcntb][0] = Y[i][t];
	    zmatb[tcntb][1] = Y[i][t+1];
	    tcntb++;
	  }
	  else {
	    zmata[tcnta][0] = Y[i][t];
	    zmata[tcnta][1] = Y[i][t+1];
	    tcnta++;
	  }
	}
	else {
	  zmatb[tcntb][0] = Y[i][t];
	  zmatb[tcntb][1] = Y[i][t+1];
	  tcntb++;
	  zmata[tcnta][0] = Y[i][t];
	  zmata[tcnta][1] = Y[i][t+1];
	  tcnta++;
	}
      }
    }

    nptra = &na;
    nptrb = &nb; 

    if ( tcnta > 0 ) mkPaTablex(zmata,tcnta,2, nptra, ((Xptr->a)[t]) );
    else na = 0;

    if ( tcntb > 0 ) mkPaTablex(zmatb,tcntb,2, nptrb, ((Xptr->b)[t]) );
    else nb = 0;

    Xptr->Na[t] = na;
    Xptr->Nb[t] = nb;

    k = 0;
    for ( i = 0; i < nb; i++ ) {
      while( ((Xptr->b)[t])[i][0] != SEnv.uvalues[k] ) k++;
      ((Xptr->Posb)[t])[i] = k;
    }
  }
  return 0;
}

// This builds the a table for a P object -- caller passes space for output table
// uses SEnv.scratch
int mkPaTablex( double **x, int n, int m, int *Tabrows , double **Table) {
  int i,j;
  int *nutableptr, nutable;
  double **utable;

  double u;
  int cnt;
 
  nutableptr = &nutable;
  utable = (SEnv.Pscrsch)->utable;
  uTabd3x( x, n, m, nutableptr, utable);

  u   = 0;
  cnt = 0;

  for ( i = 0; i < nutable; i++ ) {
    if ( i == 0 ) {
      cnt = 0;
      u = utable[i][0];
      for ( j = 0; j < 3; j++ ) Table[cnt][j] = 0;
      Table[cnt][0] = utable[i][0];
      if (  isnan( utable[i][1]) ) Table[cnt][1] = utable[i][m];
      if ( !isnan( utable[i][1]) ) Table[cnt][2] = utable[i][m];
    }
    else {
      if ( u != utable[i][0] ) {
	cnt++;
	u = utable[i][0];
	for ( j = 0; j < 3; j++ ) Table[cnt][j] = 0;
	Table[cnt][0] = utable[i][0];
	if (  isnan(utable[i][1]) ) Table[cnt][1] = utable[i][m];
	if ( !isnan(utable[i][1]) ) Table[cnt][2] = utable[i][m];
      }
      else {
	if (  isnan(utable[i][1]) ) Table[cnt][1] = Table[cnt][1] + utable[i][m];
	if ( !isnan(utable[i][1]) ) Table[cnt][2] = Table[cnt][2] + utable[i][m];
      }
    }
  }
  cnt++;

  *Tabrows = cnt;
  return 0;
}


// update an existing Pstruct using the smoothing parameter sigma
int updateP( Pstruct *Xptr, double sigma)
{
  int i,k,t;
  int type;
  int nt, na, nb;
  double X, Y;
  int ndrop, nstay;
  double num, den;
  double tmp, tmp2, sigma2;
  double wt, Dden1, Dden2, Dnum1, Dnum2;

  type         = Xptr->Type;
  nt           = Xptr->NT;
  sigma2       = pow(sigma,2);
 
  for ( t=0; t < (nt-1); t++ ) {
    na  = (Xptr->Na)[t];
    nb  = (Xptr->Nb)[t];

    for ( k = 0; k < nb; k++ ) {
      Y = ((Xptr->b)[t])[k][0];

      ((Xptr->P )[t])[k] = 0.0;
      if ( type == 0 ) {
	((Xptr->D1)[t])[k] = 0.0;
	((Xptr->D2)[t])[k] = 0.0;
      }

      num       = 0.0;
      den       = 0.0;
      Dnum1     = 0.0;
      Dden1     = 0.0;
      Dnum2     = 0.0;
      Dden2     = 0.0;

      for ( i=0; i < na; i++ ) {
	X = ((Xptr->a)[t])[i][0];

	ndrop = ((Xptr->a)[t])[i][1];
	nstay = ((Xptr->a)[t])[i][2];

	if ( (nstay+ndrop) > 0 ) {
	  tmp   = (X-Y)/sigma;
	  tmp2  = pow(tmp,2);
    	  if ( tmp2 > (2.0*720.0) ) tmp2 = 2*720.0;
	  wt    = (nstay + ndrop) * exp( -0.5 * tmp2 ); 
	  den   = den + wt;
	  if ( type == 0 ) {
	    Dden1 = Dden1 + wt * tmp2 / sigma;
	    Dden2 = Dden2 + wt * tmp2 * ( tmp2 - 3 ) / sigma2; 
	  }
	  if ( nstay > 0 ) {
	    wt    = nstay * exp( -0.5 * tmp2 ); 
	    num   = num + wt;
	    if ( type == 0 ) {
	      Dnum1 = Dnum1 + wt * tmp2 / sigma;
	      Dnum2 = Dnum2 + wt * tmp2 * ( tmp2 - 3 ) / sigma2; 
	    }
	  }
	}
      }
      if ( den > 0.0 ) {
	((Xptr->P )[t])[k] = num / den; 
	if ( type == 0 ) {
	  if ( den > 0.0 ) {
	    ((Xptr->D1)[t])[k] = ( Dnum1 - num * Dden1 / den ) / den;
	    ((Xptr->D2)[t])[k] = ( Dnum2 - 2.0 * Dnum1 * Dden1 / den - num * Dden2 / den + 2.0 * (num/den) * pow( Dden1, 2) / den) / den;
	  }
	  else {
	    ((Xptr->D1)[t])[k] = 0.0;
	    ((Xptr->D2)[t])[k] = 0.0;
	  }
	}
      }
    }
  }

  return 0;
}
