// ----------------------------------------------------------------------------
// define the cumulative beta function (used for tilting) in samon.
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "beta_cdf.h"

double beta_cdf( double x, double a, double b) {
  double tmp, betai;
       if ( x <= 0.0 ) return 0.0;
  else if ( x >= 1.0 ) return 1.0;
  if ( x == 0.0 || x == 1.0 ) tmp = 0;
  else tmp = exp(gamma_ln(a+b) - gamma_ln(a) - gamma_ln(b) + a * log(x) + b * log( 1.0 - x ));

  if ( x < ((a+1)/(a+b+2)) ) betai = tmp * beta_confrac( x, a, b ) / a;
  else betai = 1.0 - tmp * beta_confrac( 1.0-x, b, a ) / b;

  return betai;
}

double beta_confrac(double x, double a, double b) {
  int iter;
  static double eps = 1E-10;
  static int maxiter = 100;

  double am, bm, az, bz;
  double z, d, ap, bp, app, bpp, azkp;
  double z2, ap1, apb;

  am = bm = az = 1.0;

  ap1 = a + 1.0;
  apb = a + b;
  bz  = 1.0 - ( apb * x / ap1 );

  iter = 1;
  do {
    z      = iter++;
    z2     = 2.0*z;

    d      = z*(b-z)*x/(( a - 1 + z2)*(a+z2));
    ap     = az + d * am;
    bp     = bz + d * bm;

    d      = -(a+z)*(apb+z)*x/((a+z2)*(ap1+z2));
    app    = ap + d*az;
    bpp    = bp + d*bz;

    am     = ap / bpp;
    bm     = bp / bpp;

    azkp   = az;
    az     = app / bpp;
    bz     = 1;

  } while ( eps*fabs(az) <= fabs(az-azkp) && iter < maxiter ); 

  return az;
}

double gamma_ln(double x) {
  static double gammaCoef[6] = {   7.618009173E+01,   -8.650532033E+01,   2.401409822E+01,
  				  -1.231739516E+00,    1.20858003E-03,   -5.36382E-06 };
  int i;
  double tmp, series, gammaln;
  x = x - 1.0;
  tmp = x + 5.5;
  tmp = (x+0.5)*log(tmp)-tmp;
  series = 1.0;
  for ( i = 0; i < 6; i++ ) series = series + (gammaCoef[i] / ( ++x ));
  gammaln = tmp + log(2.50662827465 * series);
  return gammaln;
}
