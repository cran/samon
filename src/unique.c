#include "unique.h"
#include <stdlib.h>

double **uniqueVal( double **mat, int n, int m, int *count )
{
  int i, j;
  int cntr;
  int counter;
  double **tmpmat, **outmat;
  double *vec;

  // copy mat to a vector -- remove nan's
  vec = malloc( n * m * sizeof( double ) );
  cntr = 0;
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < m; j++ ) {
      if ( !isnan(mat[i][j]) )  vec[cntr++] = mat[i][j];
    }
  }

  // initial vector --> temporary matrix
  tmpmat = mkMatd( cntr, 2 );
  uTabxna( vec, cntr, &counter, tmpmat);

  outmat = mkMatd( counter, 2 );
  for ( i = 0; i < counter; i++ ) {
    for ( j = 0; j < 2; j++ ) {
      outmat[i][j] = tmpmat[i][j];
    }
  }
  *count = counter;

  freeMatd(tmpmat);
  free(vec);

  return outmat;
}
