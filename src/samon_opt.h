// ----------------------------------------------------------------------------
// Header for optimization of the P and Q smoothing parameters
// ----------------------------------------------------------------------------
#ifndef samon_opt_h_
#define samon_opt_h_

// optimize the smoothing parameter
int Popt(double **y, int *iter, double *optx, double *optfn);
// newton rapson method
int Pmin(int NParts, int *iter, double *optx, double *optfn);
// loss function to minimize
double lossP( double sigma, int NParts, double *deriv1, double *deriv2 );


// optimize the smoothing parameter
int Qopt(double **y, int *iter, double *optx, double *optfn);
// newton rapson method
int Qmin(int NParts, int *iter, double *optx, double *optfn);
// loss function to minimize
double lossQ(double sigma, int NParts, double *dlossqptr, double *d2lossqptr);
#endif
