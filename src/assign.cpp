#include <arprec/mp_real.h>
#include "small_inline.h"

mp_real& mp_real::operator=(const mp_real_temp& jb) {
  if (!mpr) {
    mpr = jb.mpr;
    return *this;
  }

  double n1 = mpr[0];
  double n2 = jb.mpr[0];
  //Must check to make sure that precision limits
  // for the assigned variable will not change.
  if (n1 == n2) {
    if(alloc) {
      delete [] mpr;
      mpr = jb.mpr;
    } else {
      memcpy(mpr, jb.mpr, sizeof(double) * static_cast<int>(n1));
      delete [] jb.mpr;
      /*
      for (int i = 0; i < n1; i++)
        mpr[i] = jb.mpr[i];
      */
    }
  } else if (n1 < n2) {
    // May Need to truncate some words. use rounding.
    double nw = mpr[0];
    double sgn = sign(1.0, jb.mpr[1]);
    mpr[1] = sgn * std::min(fabs(jb.mpr[1]), nw-FST_M-2.0);
    for (int i = 2; i < nw; i++)
      mpr[i] = jb.mpr[i];
    mproun(*this);

    delete [] jb.mpr;

    /*
    double tnw = this->mpr[0];
    double sgn, old_nw;
    delete [] mpr;
    mpr = jb.mpr;
    old_nw = mpr[1];
    sgn = sign(1.0, mpr[1]);
    mpr[1] = sgn * std::min(fabs(mpr[1]), double(tnw-FST_M-2.0));
    mpr[0] = tnw;
    if (old_nw == mpr[1])
      return *this;
    else {
      mproun(*this);
      return *this;
    }
    */
  } else {
    /* last case: this->mpr[0] > jb.mpr[0] 
       Need to copy to maintain same maximum number of words. */
    int i, end = int(fabs(jb.mpr[1]))+FST_M;
    for(i = 1; i < end; i++) mpr[i] = jb.mpr[i];
    delete [] jb.mpr;
    return *this;
  } 

  return *this;
}
