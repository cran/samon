// ----------------------------------------------------------------------------
// Use logistic regression to fill in intermittent data and then remove it
// ----------------------------------------------------------------------------
#ifndef _samon_FillRemove_h_
#define _samon_FillRemove_h_

int FillIn( int nv, LogisticS *logS );
void Ctime ( int t0, int *count0, int *count1, LogisticS *logS  );
int mkIM( LogisticS *logS );
#endif
