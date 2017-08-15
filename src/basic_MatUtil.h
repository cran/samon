// ----------------------------------------------------------------------------
// header for samon_MatUtil, a collection of functions to handle 1 and
// 2-dim arrays.
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// make, make and zero, zero, free a 2dim integer array
int      **mkMati(int NR, int NC);
int      **mkMatiz(int NR, int NC);
int      zeroi(int **iptr, int NR, int NC);
int      freeMati(int **iiptr);

// make, make and zero, zero, free a 2dim array of doubles
double   **mkMatd(int NR, int NC);
double   **mkMatdz(int NR, int NC);
int      zerod(double **dptr, int NR, int NC);
int      freeMatd(double **ddptr);

// copy mat N to mat M -- 2d dimensions should be appropriate
// Note this is just to make cleaner code and to hide i and j.
void cpMati( int **M, int **N, int r, int c );
void cpMatd( double **M, double **N, int r, int c );
void cpMatid( int **M, double **N, int r, int c );
void cpMatdi( double **M, int **N, int r, int c );

// sort an integer vector, a 2-d integer array on first column and
// sort an integer 2-d interger array on first 2 columns
void     qTab( int *x, int n);
void     qTab2( int **x, int n, int m);
void     qTab3( int **x, int n, int m);

// tabulate an integer vector
// output space passed by caller to uTabx
int      **uTab( int *x, int n, int *ocount);
int      uTabx( int *x, int n, int *ocount, int **optr);

// tabulates an array of ints based on the first 2 columns
// output space is passed by caller to uTab3x
int      **uTab3( int **x, int n, int m, int *Nunique);
int      uTab3x( int **x, int n, int m, int *Nunique, int **optr);

//* for doubles
// sort a vector of doubles, and sort a 2-d array on first column  
void     qTabd( double *x, int n);
void     qTabd2( double **x, int n, int m);


int uTabd3x( double **x, int n, int m, int *Nunique, double **optr);
void qTabd3( double **x, int n, int m);
double lossP( double sigma, int NParts, double *deriv1, double *deriv2 );

void qTabdna( double *x, int n);
int uTabxna( double *x, int n, int *ocount, double **optr);

int lt( double x, double y );
int eq( double x, double y );


