/*
 * src/mpreal4.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Additional routines, mostly those for which speed is less important.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

static double arprec_drand48() {
#ifdef HAVE_DRAND48
  return drand48();
#else
  /* We don't have drand48.  Use rand() to get the bits.  We call 
     rand() three times since RAND_MAX it at least 16 bits. */
  double f = 1.0 / (RAND_MAX + 1.0);
  double x = rand();
  x = x * f + rand();
  x = x * f + rand();
  return x * f;
#endif
}

void mp_real::mprand(mp_real& a)
{
  // This returns a pseudo-random MP number A between 0 and 1. 
  // Debug output starts with debug_level == 9.

  double t1;
  int nwds;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(a);
    return;
  }
  
  nwds = std::min(prec_words, int(a[0])-5);

  a[1] = nwds;
  a[2] = -1.0;
  if(mpnbt > 48) {
    for(int i =FST_M; i<nwds+FST_M;i++) {
      t1 = arprec_drand48(); 
      t1 += 0.125 * arprec_drand48();
      t1 /= 1.125; // drand48 only has 48 bit random numbers, we need mpnbt bits
      a[i] = aint(t1 * mpbdx);
    }
  } else {
    for(int i =FST_M; i<nwds+FST_M;i++) {
      a[i] = aint(arprec_drand48() * mpbdx);
    }
  }
  a[nwds+FST_M] = a[nwds+FST_M+1] = 0.0;

  // possibly (very unlikely) there are leading or trailing zeros.
  mproun(a);
  
  if(debug_level >= 9)
    cerr << "MPRAND 0" << endl;
}
