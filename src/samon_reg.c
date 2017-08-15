#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "samon_env.h"
#include "samon_reg.h"

static const
R_CMethodDef cMethods[] = {
  {"samon_boot_jk2",  (DL_FUNC) &samon_boot_jk2,   40},
  {"samon_eval",      (DL_FUNC) &samon_eval,        8},
  {"samon_gen",       (DL_FUNC) &samon_gen,        17},
  {"samon_genIM",     (DL_FUNC) &samon_genIM,      24},
  {"samon_evalIM",    (DL_FUNC) &samon_evalIM,     14},
  {"samon_ngenIMIF",  (DL_FUNC) &samon_ngenIMIF,   46},
  {NULL,NULL,0}
};

void R_init_samon(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL );
  R_useDynamicSymbols(info, TRUE);
}
