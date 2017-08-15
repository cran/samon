// ----------------------------------------------------------------------------
// a collection of functions to handle 1 and 2-dim arrays.
// ----------------------------------------------------------------------------

#include "basic_MatUtil.h"

// Make an NP by NC matrix of integers  
int** mkMati(int NR, int NC) {
  int i;
  int **iiptr, *data;

  if ( NR * NC == 0 ) return 0;

  iiptr = malloc(NR * sizeof(int*));
  data  = malloc(NR * NC * sizeof(int));
  for ( i=0; i < NR; i++ ) {
    iiptr[i] = data + (i * NC);
  }
  return iiptr;
}

// Make an NP by NC matrix of integers and zero it
int** mkMatiz(int NR, int NC) {
  int i,j;
  int **iiptr, *data;

  if ( NR * NC == 0 ) return 0;

  iiptr = malloc(NR * sizeof(int*));
  data  = malloc(NR * NC * sizeof(int));
  for ( i=0; i < NR; i++ ) {
    iiptr[i] = data + (i * NC);
    for ( j=0; j < NC; j++ ) {
      iiptr[i][j] = 0;
    }
  }
  return iiptr;
}

// zero a matrix of integers
int zeroi(int **iptr, int NR, int NC) {
  int i,j;
  for ( i = 0; i < NR; i++ ) {
    for ( j = 0; j < NC; j++ ) {
      iptr[i][j] = 0;
    }
  }
  return 0;
}

// free a matrix of integers
int freeMati(int **iiptr) {
  if ( iiptr != NULL ) {
    if ( iiptr[0] != NULL )free(iiptr[0]);
    free(iiptr);
  }
  return 0;
}

// make an NR by NC matrix of doubles
double** mkMatd(int NR, int NC) {
  int i;
  double **ddptr, *data;

  if ( NR * NC == 0 ) return 0;

  ddptr = malloc(NR * sizeof(double*) );
  data  = malloc(NR * NC * sizeof(double));
  for ( i=0; i < NR; i++ ) {
    ddptr[i] = data + (i * NC);
  }
  return ddptr;
}

// make an NR by NC matrix of doubles and zero it
double** mkMatdz(int NR, int NC) {
  int i,j;
  double **ddptr, *data;

  if ( NR * NC == 0 ) return 0;

  ddptr = malloc(NR * sizeof(double*) );
  data  = malloc(NR * NC * sizeof(double));
  for ( i=0; i < NR; i++ ) {
    ddptr[i] = data + (i * NC);
    for ( j =0; j < NC; j++ ) {
      ddptr[i][j] = 0.0;
    }
  }
  return ddptr;
}

// zero a matrix of doubles
int zerod(double **dptr, int NR, int NC) {
  int i,j;
  for ( i = 0; i < NR; i++ ) {
    for ( j = 0; j < NC; j++ ) {
      dptr[i][j] = 0.0;
    }
  }
  return 0;
}

// free a matrix of doubles
int freeMatd(double **ddptr) {
  if ( ddptr != NULL ) {
    if ( ddptr[0] != NULL ) free(ddptr[0]);
    free(ddptr);
  }
  return 0;
}

// copy mat N to mat M -- 2d dimensions should be appropriate
// Note this is just to make cleaner code and to hide i and j.
void cpMati( int **M, int **N, int r, int c ) {
  int i,j;
  for (i = 0; i < r; i++ ) {
    for (j = 0; j < c; j++ ) M[i][j] = N[i][j];
  }
}

// copy mat N to mat M -- 2d dimensions should be appropriate
// Note this is just to make cleaner code and to hide i and j.
void cpMatd( double **M, double **N, int r, int c ) {
  int i,j;
  for (i = 0; i < r; i++ ) {
    for (j = 0; j < c; j++ ) M[i][j] = N[i][j];
  }
}

// copy mat N to mat M -- 2d dimensions should be appropriate
// Note this is just to make cleaner code and to hide i and j.
void cpMatid( int **M, double **N, int r, int c ) {
  int i,j;
  for (i = 0; i < r; i++ ) {
    for (j = 0; j < c; j++ ) M[i][j] = (int) N[i][j];
  }
}

// copy mat N to mat M -- 2d dimensions should be appropriate
// Note this is just to make cleaner code and to hide i and j.
void cpMatdi( double **M, int **N, int r, int c ) {
  int i,j;
  for (i = 0; i < r; i++ ) {
    for (j = 0; j < c; j++ ) M[i][j] = (double) N[i][j];
  }
}

// sort an integer vector
void qTab( int *x, int n) {
  int i, j, p, t;

  if ( n < 2 ) return;

  p = x[ n/2 ];
  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( x[i] < p ) i++;
    while ( p < x[j] ) j--;
    if ( i >= j ) break;

    t    = x[i];
    x[i] = x[j];
    x[j] = t;
  }
  qTab( x, i );
  qTab( x + i, n - i);
}

// sort an integer matrix  on first column
void qTab2( int **x, int n, int m) {
  int i, j, p, t, k;

  if ( n < 2 ) return;

  p = x[ n/2 ][0];
  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( x[i][0] < p ) i++;
    while ( p < x[j][0] ) j--;
    if ( i >= j ) break;

    for (k=0; k<m; k++) {
      t    = x[i][k];
      x[i][k] = x[j][k];
      x[j][k] = t;
    }
  }
  qTab2( x, i, m );
  qTab2( x + i, n - i, m);
}

// sort an integer matrix on first two columns (m>1)
void qTab3( int **x, int n, int m) {
  int i, j, t, k;
  int p1, p2;

  if ( n < 2 ) return;

  p1 = x[ n/2 ][0];
  p2 = x[ n/2 ][1];

  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( x[i][0] < p1 || ( x[i][0] == p1 && x[i][1] < p2 )) i++;
    while ( p1 < x[j][0] || ( p1 == x[j][0] && p2 < x[j][1] )) j--;

    if ( i >= j ) break;

    for (k=0; k<m; k++) {
      t       = x[i][k];
      x[i][k] = x[j][k];
      x[j][k] = t;
    }
  }
  qTab3( x, i, m );
  qTab3( x + i, n - i, m);
}

// sort a vector of doubles 
void qTabd( double *x, int n) {
  int i, j;
  double p, t;

  if ( n < 2 ) return;

  p = x[ n/2 ];
  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( x[i] < p ) i++;
    while ( p < x[j] ) j--;
    if ( i >= j ) break;

    t    = x[i];
    x[i] = x[j];
    x[j] = t;
  }
  qTabd( x, i );
  qTabd( x + i, n - i);
}

// sort a array of doubles on first column
void qTabd2( double **x, int n, int m) {
  int i, j, k;
  double p, t;

  if ( n < 2 ) return;

  p = x[ n/2 ][0];
  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( x[i][0] < p ) i++;
    while ( p < x[j][0] ) j--;
    if ( i >= j ) break;

    for (k=0; k<m; k++) {
      t       = x[i][k];
      x[i][k] = x[j][k];
      x[j][k] = t;
    }
  }
  qTabd2( x, i, m );
  qTabd2( x + i, n - i, m);
}

// tabulates a vector of ints -- returns N(unique) by 2
int **uTab( int *x, int n, int *ocount) {
  int i, count;
  int **optr;

  qTab( x, n );

  // count unique values
  count = 0;
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) count = 0;
    else if ( x[i] != x[i-1] ) count++;
  }
  count++;

  optr = malloc( count * sizeof(int*));
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      count = 0;
      optr[count] = malloc( 2 * sizeof(int));
      optr[count][0] = x[i];
      optr[count][1] = 1;
    }
    else {
      if ( x[i] != x[i-1] ) {
	count++;
	optr[count] = malloc( 2 * sizeof(int));
	optr[count][0] = x[i];
	optr[count][1] = 1;
      }
      else {
	optr[count][1] = optr[count][1] + 1;
      }
    }
  }
  count++;
  *ocount = count;
  return optr;
}


// tabulates a vector of ints -- returns N(unique) by 2
// caller passes space for output
int uTabx( int *x, int n, int *ocount, int **optr) {
  int i, count;

  qTab( x, n );

  count = 0;
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      count = 0;
      optr[count][0] = x[i];
      optr[count][1] = 1;
    }
    else {
      if ( x[i] != x[i-1] ) {
	count++;
	optr[count][0] = x[i];
	optr[count][1] = 1;
      }
      else {
	optr[count][1] = optr[count][1] + 1;
      }
    }
  }
  count++;
  *ocount = count;
  return 0;
}

// tabulates an array of ints based on the first -- 2 coulumns of input array
// returns N(unique) by m+1
int **uTab3( int **x, int n, int m, int *Nunique) {

  int i, k;
  int ucnt;
  int u0, u1;
  int **optr;

  qTab3( x, n, m );
  ucnt = 0;
  u0   = 0;
  u1   = 0;
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      ucnt = 0;
      u0   = x[i][0];
      u1   = x[i][1];
    }
    else {
      if ( u0 != x[i][0] || u1 != x[i][1] ) {
	ucnt++;
	u0 = x[i][0];
	u1 = x[i][1];
      }
    }
  } 
  ucnt++;
  optr = malloc( ucnt * sizeof( int* ) );
  for ( i = 0; i < ucnt; i++ ) {
    optr[i] = malloc( (m+1) * sizeof( int ));
  } 
  //  and again to populate the new array
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      ucnt = 0;
      u0   = x[i][0];
      u1   = x[i][1];
      for ( k = 0; k < m; k++ ) {
	optr[ucnt][k] = x[i][k];
      }
      optr[ucnt][m] = 1;
    }
    else {
      if ( u0 != x[i][0] || u1 != x[i][1] ) {
	ucnt++;
	u0 = x[i][0];
	u1 = x[i][1];
	for ( k = 0; k < m; k++ ) {
	  optr[ucnt][k] = x[i][k];
	}
	optr[ucnt][m] = 1;
      }
      else {
	optr[ucnt][m]++;
      }
    }
  } 
  ucnt++;
  *Nunique = ucnt;
  return optr;
}

// tabulates an array of ints based on the first 2 coulumns of input array
// returns N(unique) by m+1 -- caller passes space for output
int uTab3x( int **x, int n, int m, int *Nunique, int **optr) {

  int i,k;
  int ucnt;
  int u0, u1;

  qTab3( x, n, m );

  ucnt = 0; 
  u0   = 0;
  u1   = 0;
  //  and again to populate the new array
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      ucnt = 0;
      u0   = x[i][0];
      u1   = x[i][1];
      for ( k = 0; k < m; k++ ) {
	optr[ucnt][k] = x[i][k];
      }
      optr[ucnt][m] = 1;
    }
    else {
      if ( u0 != x[i][0] || u1 != x[i][1] ) {
	ucnt++;
	u0 = x[i][0];
	u1 = x[i][1];
	for ( k = 0; k < m; k++ ) {
	  optr[ucnt][k] = x[i][k];
	}
	optr[ucnt][m] = 1;
      }
      else {
	optr[ucnt][m]++;
      }
    }
  } 
  ucnt++;
  *Nunique = ucnt;
  return 0;
}

int lt( double x, double y )
{
  if (isnan(y)) return !isnan(x);
  return x < y;
}

int eq( double x, double y )
{
  if (isnan(y)) return isnan(x);
  return x == y;
}


// sort a double matrix on first two columns (m>0)
void qTabd3( double **x, int n, int m) {
  int i, j, k;
  double t;
  double p1, p2;

  if ( n < 2 ) return;

  p1 = x[ n/2 ][0];
  p2 = x[ n/2 ][1];

  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( lt(x[i][0],p1) || ( eq(x[i][0],p1) && lt(x[i][1],p2) )) i++;
    while ( lt(p1,x[j][0]) || ( eq(p1,x[j][0]) && lt(p2,x[j][1]) )) j--;

    if ( i >= j ) break;

    for (k=0; k<m; k++) {
      t       = x[i][k];
      x[i][k] = x[j][k];
      x[j][k] = t;
    }
  }
  qTabd3( x, i, m );
  qTabd3( x + i, n - i, m);
}

// tabulates an array of doubles based on the first 2 coulumns of the input array
// returns N(unique) by m+1 -- caller passes space for output
int uTabd3x( double **x, int n, int m, int *Nunique, double **optr) {

  int i,k;
  int ucnt;
  double u0, u1;

  qTabd3( x, n, m );

  ucnt = 0;
  u0   = 0;
  u1   = 0; 
  //  and again to populate the new array
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      ucnt = 0;
      u0   = x[i][0];
      u1   = x[i][1];
      for ( k = 0; k < m; k++ ) {
	optr[ucnt][k] = x[i][k];
      }
      optr[ucnt][m] = 1;
    }
    else {
      if ( !eq(u0,x[i][0]) || !eq(u1,x[i][1]) ) {
	ucnt++;
	u0 = x[i][0];
	u1 = x[i][1];
	for ( k = 0; k < m; k++ ) {
	  optr[ucnt][k] = x[i][k];
	}
	optr[ucnt][m] = 1;
      }
      else {
	optr[ucnt][m]++;
      }
    }
  } 
  ucnt++;
  *Nunique = ucnt;
  return 0;
}

// -----------------------------------------------------

// sort a vector of doubles 
void qTabdna( double *x, int n) {
  int i, j;
  double p, t;

  if ( n < 2 ) return;

  p = x[ n/2 ];
  for ( i=0, j=n-1 ;; i++, j-- ) {
    while ( lt(x[i],p) ) i++;
    while ( lt(p,x[j]) ) j--;
    if ( i >= j ) break;

    t    = x[i];
    x[i] = x[j];
    x[j] = t;
  }
  qTabdna( x, i );
  qTabdna( x + i, n - i);
}


// tabulates a vector of doubles -- returns N(unique) by 2
// caller passes space for output
int uTabxna( double *x, int n, int *ocount, double **optr) {
  int i, count;

  qTabdna( x, n );

  count = 0; 
  for ( i = 0; i < n; i++ ) {
    if ( i == 0 ) {
      count = 0;
      optr[count][0] = x[i];
      optr[count][1] = 1;
    }
    else {
      if ( eq(x[i],x[i-1]) ) {
	optr[count][1] = optr[count][1] + 1;
      }
      else {
	count++;
	optr[count][0] = x[i];
	optr[count][1] = 1;
      }
    }
  }
  count++;
  *ocount = count;
  return 0;
}
// -----------------------------------------------------
