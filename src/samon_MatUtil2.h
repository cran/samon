// ----------------------------------------------------------------------------
// Some matrix functions to support logistic regression
// Inverts positive definite matrix. 
// ----------------------------------------------------------------------------
#ifndef _samon_MatUtil2_h_
#define _samon_MatUtil2_h_
int cholesky( double **M, int n, double **L, int ** Tmodel);
int invLTri( double **L, int n, double **Linv);
int MatMult( double **M1, int nr1, int nc1, double **M2, int nr2, int nc2, double **Mo );
int Matinv( double **M, int n, double **Minv, double **L, double **Linv, int *Npivits, int **Tmodel, int **diag );
#endif
