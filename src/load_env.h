// ----------------------------------------------------------------------------
// Header for loading the samon environment
// ----------------------------------------------------------------------------
#ifndef load_env_h_
#define load_env_h_
// load the environment
int load_env ( 
   int            N0,   int           NT,  
   int         seed0,   double    startp,   double     HSigp,
   double     startq,   double     HSigq,
   double         lb,   double        ub,   double     zeta1,   double     zeta2,
   int        NParts,
   int      NSamples,
   int       MaxIter,   double   FAconvg,   double   FRconvg,   double  SAconvg,
   int       nunique,   double  **uvalues,  double    smallv 
	       );

// reload when N0 is decreased for example
int reload_env ( int N0, double startp, double startq );

int free_env(void);

// needed for creating the right environment
int partition( int ***Part, int NParts, int N0 );
#endif
