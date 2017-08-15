// ----------------------------------------------------------------------------
// Header for IF functions
// ----------------------------------------------------------------------------

#ifndef IF_fun_beta_h_
#define IF_fun_beta_h_

// some functions to build matrices

int mkInt( Pstruct *Pptr, Qstruct *Qptr);
int  mkH(Qstruct *Qptr); 
int mkT2(void);
int mkB3(void);

int mkBMat(Qstruct *Qptr);
int mkC(void);
int mkUU(int, double, double, double, double*, double*, double*, double*, Qstruct*, int, double**);

//int mkPre(void);
int mkPre(Pstruct *Pptr, Qstruct *Qptr);

double beta_cdf( double x, double a, double b);

int mkQ0(Pstruct *X);
int mkPQMat( Pstruct *Pptr, Qstruct *Qptr, int nr, int nc, int nt, double alpha );

int update(Qstruct* Qptr, double alpha);
int updateT(Qstruct* Qptr, double alpha);
int mkV(int **Y, int nr, int nt, int n0);

int IF_fun(double **y, double in_sigmap, double in_sigmaq, int tp1, int rep, double **OOptr, int Nalpha, double *alphalist, int ifi, double **ifiptr );
#endif
