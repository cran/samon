
void samon_gen ( 
  int               *N0,   int               *NT,    double          *Mat, 

  double      *outPmatM,   double      *outQmatM,

  int            *seed0,
  double        *startp,   double         *HSigp,
  double        *startq,   double         *HSigq,
  int           *NParts,
  int         *NSamples,
  int          *MaxIter,
  double       *FAconvg,   double       *FRconvg,   double       *SAconvg,
  double        *Sample
		 );

void samon_eval ( 
   int            *N0,    int            *NT,  double        *Mat, 
   double     *outmat,
   int        *nsigma, double    *SigmaList,
   int        *NParts,    int          *PQtp
		  );

void samon_boot_jk2 ( 
  int               *N0,   int               *NT,    double          *Mat, 

  double      *outPmatM,   double      *outQmatM,   double     *alphamatM,
  double    *outPmatMjk,   double    *outQmatMjk,   double   *alphamatMjk,
  double      *outPmatS,   double      *outQmatS,   double     *alphamatS,
  double    *outPmatSjk,   double    *outQmatSjk,   double   *alphamatSjk,

  int            *seed0,
  double        *startp,   double         *HSigp,
  double        *startq,   double         *HSigq,
  double            *lb,   double            *ub,
  double         *zeta1,   double         *zeta2,
  int           *NParts,
  int         *NSamples,
  int          *MaxIter,
  double       *FAconvg,   double       *FRconvg,   double       *SAconvg,
  int           *Nalpha,   double     *alphaList,
  int              *Mjk,   int              *Sjk,
  int          *retIFiM,   int          *retIFiS,
  double       *IFvalsM,   double       *IFvalsS,
  int          *lenTfun,   double          *Tfun
		      );

void samon_genIM ( 
   int            *N0,  int           *NT,  double       *Mat, 
   int        *nmodel,  int      *inmodel,
   double       *FMat,  double    *LEstsM,  int      *ModelsM,
   double   *outPmatM,  double  *outQmatM,
   int         *seed0,  int        *seed1,
   double     *startp,  double     *HSigp,
   double     *startq,  double     *HSigq,
   int        *NParts,
   int      *NSamples, 
   int       *MaxIter,  double   *FAconvg,  double   *FRconvg,  double  *SAconvg,
   int        *NFills,
   double     *Sample
		   );

void samon_evalIM ( 
   int            *N0,  int            *NT,  double        *Mat, 
   int        *nmodel,  int       *inmodel,
   double       *FMat,  double     *LEstsM,  int      *ModelsM,
   double     *outmat,
   int        *nsigma,  double  *SigmaList,
   int         *seed1,
   int        *NParts,
   int          *PQtp
		    );

void samon_ngenIMIF ( 
   int            *N0,  int           *NT,  double       *Mat, 
   int        *nmodel,  int      *inmodel,
   double      *FMatM,  double    *LEstsM,  int      *ModelsM,
   double      *FMatS,  double    *LEstsS,  int      *ModelsS,
   double   *outPmatM,  double  *outQmatM,  double *alphamatM,
   double   *outPmatS,  double  *outQmatS,  double *alphamatS,
   int         *seed0,  int        *seed1,
   double     *startp,  double     *HSigp,
   double     *startq,  double     *HSigq,
   double         *lb,  double        *ub,
   double      *zeta1,  double     *zeta2,
   int        *NParts,
   int      *NSamples, 
   int       *MaxIter,  double   *FAconvg,  double   *FRconvg,  double  *SAconvg,
   int        *Nalpha,  double *alphaList,
   int        *NFills,
   int       *retIFiM,  int      *retIFiS,
   double    *IFvalsM,  double   *IFvalsS,
   int     *retSample,
   double     *Sample,
   int      *retFMatM,  int     *retFMatS, 
   int       *lenTfun,  double      *Tfun
		      );

