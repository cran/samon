// ----------------------------------------------------------------------------
// Header for generation functions
// ----------------------------------------------------------------------------

#ifndef Gen_fun_h_
#define Gen_fun_h_

// #include <R.h>

// generate data from P and Q structures
int Gen_fun(double **y, int n0, int nt, double **sampy, int Nsamp, double sigp, double sigq ) ;

// initialize the random number generator from a seed
int  seed_sgen(unsigned long long int inseed);

// generate a double
double sgen(void);
#endif
