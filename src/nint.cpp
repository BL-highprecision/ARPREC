/*
 * src/mpreal.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include <arprec/mp_real.h>
#include "small_inline.h"

void mp_real::mpnint(const mp_real& a, mp_real& b, int prec_words)
{
  /**
   * Rounds to nearest integer.  If prec_words is low, 
   * so that an integer cannot be represented in 
   * fully, this may return incorrect results.
   * In other words, in large integer applications,
   * care must be given to setting prec_words sufficiently high.
   * if prec_words is too low, do not expect mpnint to produce
   * the correctly rounded integer.
   *
   * 
   */
  int na = int(std::abs(a[1]));
  int nb = std::min(int(b[0])-5, std::min(prec_words, na));
  double mag_up = 0.0;
  double exp = a[2];
  double sgn = sign(1.0, a[1]);

  if(na == 0 || nb == 0) {
    zero(b);
    return;
  }
  
  if(exp - na +1 < 0) {
    // we must truncate to exp+1 words
    int nout = std::min(nb, int(exp)+1);
    if(nout < 0) {
      // magnitude far less than 1.
      zero(b);
      return;
    }
    if(!nout) {
      // might round up to one.
      if(a[FST_M] >= mpbdx/2.0) {
        b[1] = sgn * 1.0;
        b[2] = 0.0;
        b[3] = 1.0;
        b[4] = 0.0;
        return;
      } else {
        zero(b);
        return;
      }
    }
    if(na > nout) {
      // rounding required
      if(a[FST_M + nout] >= mpbdx/2.0) {
        mag_up = 1.0;
      }
    }
    nb = nout;
  }
  b[nb+FST_M] = b[nb + FST_M+1] = 0.0;
  for(int i = nb+FST_M-1; i >= FST_M; i--) b[i] = a[i];
  b[1] = sgn * nb;
  b[2] = a[2];
  if(mag_up == 1.0) {
    // add one (or subtract one if negative).
    mp_real sk0,f;
    f[1] = 1;
    f[2] = 0;
    f[3] = 1;
    f[4] = 0;
    if(sgn > 0)
      mpadd(b, f, sk0, prec_words);
    else
      mpsub(b, f, sk0, prec_words);
    mpeq(sk0, b, prec_words);
  }
  return;
}

