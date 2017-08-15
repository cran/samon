// ----------------------------------------------------------------------------
// partition: calculate the start and stop observations that define a partiton.
//              makes a table to be put in the environment.
//   update : change interface / eliminate environment
// ----------------------------------------------------------------------------

#include "samon_env.h"
#include "load_env.h"
#include "basic_MatUtil.h"

//  ------------------------------------------------------------------------------
int partition( int ***Part, int NParts, int N0 )
{
  int i;
  int cnt, step, rem;

  *Part = mkMati( NParts, 2 );

  rem = N0 % NParts;
  cnt = 0;
  for ( i = 0; i < NParts; i++ ) {

    if ( i < (NParts - rem) ) step = floor(N0/NParts) - 1;
    else step = floor(N0/NParts);

    (*Part)[i][0] = cnt;
    (*Part)[i][1] = cnt + step;

    cnt = cnt + step + 1; 
  }
  return 0;
}
