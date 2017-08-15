// ----------------------------------------------------------------------------
// Create, populate, maintain and delete Qstruct objects
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "basic_MatUtil.h"
#include "Qstruct.h"

// -----------------------------------------------------
// -----------------------------------------------------

int posc( Qstruct *Xptr, int t, double Y )
{
  int i, nc;

  nc = (Xptr->Nc)[t];
  if ( nc == 0 ) return -1;
  if ( Y <= ((Xptr->Qc)[t])[0] ) return 0;
  for ( i = 1; i < nc; i++ ) {
    if ( (((Xptr->Qc)[t])[i-1] < Y) && (Y <= ((Xptr->Qc)[t])[i] ) ) return i;
  }
  return -1;
}

int posr( Qstruct *Xptr, int t, double Y )
{
  int i, nr;

  nr = (Xptr->Nr)[t];
  if ( nr == 0 ) return -1;
  if ( Y <= ((Xptr->Qr)[t])[0] ) return 0;
  for ( i = 1; i < nr; i++ ) {
    if ( (((Xptr->Qr)[t])[i-1] < Y) && (Y <= ((Xptr->Qr)[t])[i] ) ) return i;
  }
  return -1;
}


// -----------------------------------------------------
// distruct a Qstruct
int Qdestruct( Qstruct *X ) {
  int t, type;
  int nt;

  type = X->Type;
  nt   = X->NT;

  for ( t = 0; t < (nt-1); t++ ) {

    freeMatd( X->a[t] );
    freeMatd( X->b[t] );

      freeMatd( X->Q[t]  );
      freeMatd( X->CQ[t] );
      if ( type == 0 ) {
	freeMatd( X->D1[t]   );
	freeMatd( X->D2[t]   );
	freeMatd( X->DCQ1[t] );
	freeMatd( X->DCQ2[t] );
      }

      if ( type == 1 ) {
	freeMatd( X->TQ[t]   );
	freeMatd( X->H [t]   );
	freeMatd( X->IFB[t]  );
      }

      free( X->Qr[t]   );
      free( X->Qc[t]   );
      free( X->Posr[t] );
      free( X->Posc[t] );
      free( X->Qd[t]  );
  }
  free( X->Na   );
  free( X->Nb   );

  free( X->acnt );
  free( X->bcnt );
  free( X->Nr   );
  free( X->Nc   );

  free( X->Qr   );
  free( X->Qc   );
  free( X->Posr );
  free( X->Posc );
  free( X->Qd   );

  free( X->a    );
  free( X->b    );
  free( X->Q    );
  free( X->CQ   );
  if ( type == 0 ) {
    free( X->D1   );
    free( X->D2   );
    free( X->DCQ1 );
    free( X->DCQ2 );
  }
  if ( type == 1 ) {
    free( X->TQ   );
    free( X->H    );
    free( X->IFB  );
  }

  free( X       );

  return 0;
} 

// -----------------------------------------------------
// -----------------------------------------------------

// update an existing Qstruct using the smoothing parameter sigma
int updateQ( Qstruct *Xptr, double sigma )
{
  int i, j, k, t;
  int type;

  int pos;
  double wt;
  double A, dA, d2A;
  double B, dB, d2B;
  double sigma2;
  double z2;
  double BYd, AY0d, AY1d;
  double nnd;

  int nt, na, nr, nc;

  type    = Xptr->Type;
  nt      = Xptr->NT;
  sigma2  = pow(sigma,2.0); 

  for ( t=0; t < (nt-1); t++ ) {

    na = (Xptr->Na)[t];
    nr = (Xptr->Nr)[t];
    nc = (Xptr->Nc)[t];

    for ( i = 0; i < nr; i++ ) {

      // set to zero
      for ( j = 0; j < nc; j++ ) {
	((Xptr->Q   )[t])[i][j] = 0.0;
	((Xptr->CQ  )[t])[i][j] = 0.0;
	if ( type == 0 ) {
	  ((Xptr->D1  )[t])[i][j] = 0.0;
	  ((Xptr->D2  )[t])[i][j] = 0.0;
	  ((Xptr->DCQ1)[t])[i][j] = 0.0;
	  ((Xptr->DCQ2)[t])[i][j] = 0.0;
	}
      }

      BYd  = ((Xptr->Qr)[t])[i];

      B   = 0.0;
      dB  = 0.0;
      d2B = 0.0;
      for ( k = 0; k < na; k++ ) {
	AY0d  = ((Xptr->a)[t])[k][0];
	AY1d  = ((Xptr->a)[t])[k][1];
	nnd   = ((Xptr->a)[t])[k][2];
	pos   = posc(Xptr,t,AY1d);

	z2    = pow( ( AY0d - BYd ) / sigma, 2.0);
	if ( z2 > 2*720 ) z2 = 2*720;
	wt  = nnd * exp( -0.5 * z2 );
	((Xptr->Q )[t])[i][ pos ] = ((Xptr->Q )[t])[i][ pos ] + wt;
	if ( type == 0 ) {
	  ((Xptr->D1)[t])[i][ pos ] = ((Xptr->D1)[t])[i][ pos ] + ( wt * z2 / sigma );
	  ((Xptr->D2)[t])[i][ pos ] = ((Xptr->D2)[t])[i][ pos ] + ( wt * z2 * ( z2 - 3.0 ) / sigma2 );
	}

	B   = B + wt;
	if ( type == 0 ) {
	  dB  = dB + ( wt * z2 / sigma );
	  d2B = d2B + ( wt * z2 * ( z2 - 3.0 ) / sigma2 );
	}

      }

      for (k = 0; k < nc; k++ ) {
	if ( B > 0.0 ) {
	  A   = ((Xptr->Q )[t])[i][k];
	  if ( type == 0 ) {
	    dA  = ((Xptr->D1)[t])[i][k];
	    d2A = ((Xptr->D2)[t])[i][k];
	  }

	  if ( type == 0 ) {
	    ((Xptr->D1)[t])[i][k] = ( dA / B - (A * dB / B) / B  );
	    ((Xptr->D2)[t])[i][k] = ( d2A / B - 2.0 * (dA/B) * (dB/B) - (A/B) * (d2B/B) + 2.0 * ( A/B ) * pow( dB/B, 2.0) );
	  }
	  ((Xptr->Q )[t])[i][k] = A / B;
	}
      }

      // cumulative
      if ( type == 0 ) {
	for (k = 0; k < nc; k++ ) {
	  if ( k == 0 ) {
	    ((Xptr->CQ  )[t])[i][k] = ((Xptr->Q )[t])[i][k];
	    ((Xptr->DCQ1)[t])[i][k] = ((Xptr->D1)[t])[i][k];
	    ((Xptr->DCQ2)[t])[i][k] = ((Xptr->D2)[t])[i][k];
	  }
	  else {
	    ((Xptr->CQ  )[t])[i][k] = ((Xptr->Q )[t])[i][k] + ((Xptr->CQ  )[t])[i][k-1];
	    ((Xptr->DCQ1)[t])[i][k] = ((Xptr->D1)[t])[i][k] + ((Xptr->DCQ1)[t])[i][k-1];
	    ((Xptr->DCQ2)[t])[i][k] = ((Xptr->D2)[t])[i][k] + ((Xptr->DCQ2)[t])[i][k-1];
	  }
	}
      }
      else {
	for (k = 0; k < nc; k++ ) {
	  if ( k == 0 ) {
	    ((Xptr->CQ  )[t])[i][k] = ((Xptr->Q )[t])[i][k];
	  }
	  else {
	    ((Xptr->CQ  )[t])[i][k] = ((Xptr->Q )[t])[i][k] + ((Xptr->CQ  )[t])[i][k-1];
	  }
	}
      }
    }
  }

  return 0;
}

// -----------------------------------------------------
// -----------------------------------------------------

// Create a Qstruct but leave the population til later when the data becomes available
Qstruct *Qinit0( int N0, int NT, int size, int type ) {

  int t,i,j;
  int MX;

  Qstruct *Xptr;

  if ( 0 < size && size < N0 ) MX = size;
  else MX = N0;

  Xptr = malloc(sizeof(Qstruct));
  Xptr->NT   = NT;
  Xptr->N    = N0;
  Xptr->Type = type;

  Xptr->Na   = malloc( NT * sizeof( int     ));
  Xptr->Nb   = malloc( NT * sizeof( int     ));
  Xptr->Nr   = malloc( NT * sizeof( int     ));
  Xptr->Nc   = malloc( NT * sizeof( int     ));

  Xptr->acnt = malloc( NT * sizeof( int     ));
  Xptr->bcnt = malloc( NT * sizeof( int     ));

  Xptr->Qr   = malloc( NT * sizeof( double* ));
  Xptr->Qc   = malloc( NT * sizeof( double* ));

  Xptr->Posr = malloc( NT * sizeof( int*    ));
  Xptr->Posc = malloc( NT * sizeof( int*    ));

  Xptr->Qd   = malloc( NT * sizeof( double* ));

  for ( i=0; i < NT; i++ ) {
    Xptr->Na[i]   = 0;
    Xptr->Nb[i]   = 0;
    Xptr->Nr[i]   = 0;
    Xptr->Nc[i]   = 0;
    Xptr->acnt[i] = 0;
    Xptr->bcnt[i] = 0;
  }

  (Xptr->a   ) = malloc( NT * sizeof( int**    ));
  (Xptr->b   ) = malloc( NT * sizeof( int**    ));
  (Xptr->Q   ) = malloc( NT * sizeof( double** ));
  (Xptr->CQ  ) = malloc( NT * sizeof( double** ));

  if ( type == 0 ) {
    (Xptr->D1  ) = malloc( NT * sizeof( double** ));
    (Xptr->D2  ) = malloc( NT * sizeof( double** ));
    (Xptr->DCQ1) = malloc( NT * sizeof( double** ));
    (Xptr->DCQ2) = malloc( NT * sizeof( double** ));
  }
  if ( type == 1 ) {
    (Xptr->TQ  ) = malloc( NT * sizeof( double** ));
    (Xptr->H   ) = malloc( NT * sizeof( double** ));
    (Xptr->IFB ) = malloc( NT * sizeof( double** ));
  }

  for ( t = 0; t < NT-1; t++ ) {

    Xptr->a[t] = mkMatd( N0, 3 );
    Xptr->b[t] = mkMatd( N0, 3 );

  // -----------------------------------------------------

    (Xptr->Qr)[t] = malloc( MX * sizeof( double ));
    (Xptr->Qc)[t] = malloc( MX * sizeof( double ));
    (Xptr->Qd)[t] = malloc( MX * sizeof( double ));

    (Xptr->Posr)[t] = malloc( MX * sizeof( int ));
    (Xptr->Posc)[t] = malloc( MX * sizeof( int ));

    Xptr->Q[t]     = mkMatd( MX, MX );
    Xptr->CQ[t]    = mkMatd( MX, MX );
    if ( type == 0 ) {
      Xptr->D1[t]    = mkMatd( MX, MX );
      Xptr->D2[t]    = mkMatd( MX, MX );
      Xptr->DCQ1[t]  = mkMatd( MX, MX );
      Xptr->DCQ2[t]  = mkMatd( MX, MX );
    }
    if ( type == 1 ) {
      Xptr->TQ[t]    = mkMatd( MX, MX );
      Xptr->H[t]     = mkMatd( MX, MX );
      Xptr->IFB[t]   = mkMatd( MX, MX );
    }

    for ( i = 0; i < MX; i++ ) {
      for ( j = 0; j < MX; j++ ) {
	((Xptr->Q   )[t])[i][j] = 0.0;
	((Xptr->CQ  )[t])[i][j] = 0.0;
	if ( type == 0 ) {
	  ((Xptr->D1  )[t])[i][j] = 0.0;
	  ((Xptr->D2  )[t])[i][j] = 0.0;
	  ((Xptr->DCQ1)[t])[i][j] = 0.0;
	  ((Xptr->DCQ2)[t])[i][j] = 0.0;
	}
	if ( type == 1 ) {
	  ((Xptr->TQ  )[t])[i][j] = 0.0;
	  ((Xptr->H   )[t])[i][j] = 0.0;
	  ((Xptr->IFB )[t])[i][j] = 0.0;
	}
      }
    }
  }
  return Xptr;
}

// -----------------------------------------------------
// -----------------------------------------------------

// Populate an existing Qstruct -- i.e. refill the a and b tables 
int Qinit1( Qstruct *Xptr, double **Y, int N0, int NT, int start_cut, int stop_cut, int type ) {

  int t,i,j, k;
  //  int rc;
  int tcnta, tcntb, tcntp0, tcntp1;
  double **zmata, **zmatb, *zmatp0, *zmatp1;
  int na, nb, np0, np1;
  int *nptra, *nptrb, *nptrp0, *nptrp1;
  double **tablep0, **tablep1;

  int nr, nc;

  for ( t=0; t < NT; t++ ) {
    Xptr->Na[t]   = 0;
    Xptr->Nb[t]   = 0;
    Xptr->Nr[t]   = 0;
    Xptr->Nc[t]   = 0;
    Xptr->acnt[t] = 0;
    Xptr->bcnt[t] = 0;
  }

  zmata   = (SEnv.Qscrsch)->zmata;
  zmatb   = (SEnv.Qscrsch)->zmatb;
  zmatp0  = (SEnv.Qscrsch)->zmatp0;
  zmatp1  = (SEnv.Qscrsch)->zmatp1;

  tablep0 = (SEnv.Qscrsch)->tablep0;
  tablep1 = (SEnv.Qscrsch)->tablep1;

  for ( t = 0; t < NT-1; t++ ) {

    tcnta  = 0;
    tcntb  = 0;
    tcntp0 = 0;   // 1-d time-point t    b only
    tcntp1 = 0;   // 1-d time-point t+1  both a and b
    for ( i = 0; i < N0; i++ ) {
      if ( !isnan(Y[i][t]) ) {
	if ( type == 0 ) {
	  if ( !isnan(Y[i][t+1]) ) { 
	    if ( start_cut <= i && i <= stop_cut ) {
	      zmatb[tcntb][0] = Y[i][t];
	      zmatb[tcntb][1] = Y[i][t+1];
	      tcntb++;
	      zmatp0[tcntp0]  = Y[i][t];     // bs only
	      tcntp0++;
	      zmatp1[tcntp1]  = Y[i][t+1];
	      tcntp1++;
	    }
	    else {
	      zmata[tcnta][0] = Y[i][t];
	      zmata[tcnta][1] = Y[i][t+1];
	      tcnta++;
	      zmatp1[tcntp1]  = Y[i][t+1];
	      tcntp1++;
	    }
	  }
	}
	else {
	  zmatb[tcntb][0] = Y[i][t];
	  zmatb[tcntb][1] = Y[i][t+1];   // might be -1 ok
	  tcntb++;
	  zmatp0[tcntp0]  = Y[i][t];     // bs only
	  tcntp0++;

	  if ( !isnan(Y[i][t+1]) ) {
	    zmatp1[tcntp1]  = Y[i][t+1];
	    tcntp1++;

	    zmata[tcnta][0] = Y[i][t];
	    zmata[tcnta][1] = Y[i][t+1];
	    tcnta++;
	  }
	}
      }
    }

    nptra  = &na;
    nptrb  = &nb; 
    nptrp1 = &np1;
    nptrp0 = &np0;

    if (  tcnta  > 0 ) uTabd3x(zmata,  tcnta,  2,   nptra, ((Xptr->a)[t]) );
    else na = 0;

    if (  tcntb  > 0 ) uTabd3x(zmatb,  tcntb,  2,   nptrb, ((Xptr->b)[t]) );
    else nb = 0;

    if (  tcntp0 > 0 ) uTabxna (zmatp0, tcntp0,     nptrp0, tablep0);
    else np0 = 0;

    if (  tcntp1 > 0 ) uTabxna (zmatp1, tcntp1,     nptrp1, tablep1);
    else np1 = 0;

    Xptr->Na[t] = na;
    Xptr->Nb[t] = nb;

    Xptr->acnt[t] = tcnta;
    Xptr->bcnt[t] = tcntb;

  // -----------------------------------------------------

    nr = np0;
    nc = np1;
    (Xptr->Nr)[t] = nr;
    (Xptr->Nc)[t] = nc;

    k = 0;
    for ( i = 0; i < nr; i++ ) {
      ((Xptr->Qr)[t])[i] = tablep0[i][0];
      while( tablep0[i][0] != SEnv.uvalues[k] ) k++;
      ((Xptr->Posr)[t])[i] = k;
    }
    k = 0;
    for ( i = 0; i < nc; i++ ) {
      ((Xptr->Qc)[t])[i] = tablep1[i][0];
      while( tablep1[i][0] != SEnv.uvalues[k] ) k++;
      ((Xptr->Posc)[t])[i] = k;
      ((Xptr->Qd)[t])[i] = ( ( double ) tablep1[i][1] ) / ( ( double ) tcntp1 );
    }

    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < nc; j++ ) {
	((Xptr->Q   )[t])[i][j] = 0.0;
	((Xptr->CQ  )[t])[i][j] = 0.0;
	if ( type == 0 ) {
	  ((Xptr->D1  )[t])[i][j] = 0.0;
	  ((Xptr->D2  )[t])[i][j] = 0.0;
	  ((Xptr->DCQ1)[t])[i][j] = 0.0;
	  ((Xptr->DCQ2)[t])[i][j] = 0.0;
	}
	if ( type == 1 ) {
	  ((Xptr->TQ  )[t])[i][j] = 0.0;
	  ((Xptr->H   )[t])[i][j] = 0.0;
	  ((Xptr->IFB )[t])[i][j] = 0.0;
	}
      }
    }
  }
return 0;
}
