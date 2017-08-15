#include "samon_env.h"
#include "basic_MatUtil.h"
#include "samon_MatUtil2.h"

// Invert a lower triangle matrix in L: result in Linv
int invLTri( double **L, int n, double **Linv) {
  int i,j,row;
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      if ( i == j ) Linv[i][j] = 1.0;
      else Linv[i][j] = 0.0;
    }
  }
  for ( row = 0; row < n; row++ ) {
    for ( j = 0; j < (row+1); j++ ) Linv[row][j] = Linv[row][j] / L[row][row]; 
    for ( i = row+1; i < n; i++ ) {
      for ( j = 0; j < row+1; j++) {
	Linv[i][j] = Linv[i][j] - Linv[row][j] * L[i][row];
      }
    }
  }
  return 0;
}

// cholesky lower triangular matrix of a pos-def matrix M 
int cholesky( double **M, int n, double **L, int **Tmodel ) {
  int i,j,k;
  double tmp;
  double leps = 10E-10;

  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      L[i][j] = 0.0;
    }
  }
  
  for ( i = 0; i < n; i++ ) {
    if ( Tmodel[i][0] == 1 ) {
      for ( j = 0; j < (i+1); j++ ) {
	if ( Tmodel[j][0] == 1 ) {
	  tmp = 0.0;
	  for ( k = 0; k < j; k++ ) {
	    tmp = tmp + L[i][k] * L[j][k];
	  }
	  if ( i == j ) {
	    if ( (M[i][i] - tmp) < 0.0 ) L[i][j] = 0.0;
	    else L[i][j] = sqrt(M[i][i] - tmp);
	  }
	  else {
	    if ( fabs(L[j][j]) < leps ) {
	    }
	    else {
	      L[i][j] = (1.0 / ( L[j][j]) * (M[i][j] - tmp));
	    }
	  }
	}
      }
    }
    else {
      L[i][i] = 1;
    }
  }
  return 0;
}

// inverse of a positive def matrix M
int Matinv( double **M, int n, double **Minv, double **L, double **Linv, int *Npivits, int **Tmodel, int **diag ) {
  int i,j, k;
  static double leps = 10E-14;

  cholesky( M, n, L, Tmodel);
  for ( i = 0; i < n; i++ )  {
    if ( Tmodel[i][0] == 0 ) {
      for ( j = 0; j < n; j++ ) {
	L[i][j] = 0.0;
	L[j][i] = 0.0;
      }
      L[i][i] = 1.0;
    }
  }

  // now check the diagonal
  *Npivits = n;
  for ( i = 0; i < n;  i++ ) {
    diag[i][0] = 0;
    if ( fabs(L[i][i]) < leps || isnan(L[i][i]) || isinf(L[i][i]) ) {
      *Npivits = *Npivits - 1;
           if ( fabs(L[i][i] ) < leps ) diag[i][0] = 1;
      else if ( isnan(L[i][i])        ) diag[i][0] = 2;
      else if ( isinf(L[i][i])        ) diag[i][0] = 3;
    }
  }

  if ( *Npivits < n ) return 1;

  // invert the L
  invLTri(L,n,Linv);

  // multiply to get inverse of M
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ) {
      Minv[i][j] = 0.0;
      for ( k = 0; k < n; k++ ) { 
	Minv[i][j] = Minv[i][j] + Linv[k][i] * Linv[k][j];
      }
    }
  } 
  return 0;
}

int MatMult( double **M1, int nr1, int nc1, double **M2, int nr2, int nc2, double **Mo )
{
  int i,j,k;
  for ( i = 0; i < nr1; i++ ) {
    for ( j = 0; j < nc2; j++ ) {
      Mo[i][j] = 0.0;
      for ( k = 0; k < nc1; k++ ) {
	Mo[i][j] = Mo[i][j] + M1[i][k] * M2[k][j];
      }
    }
  }
  return 0;
}
