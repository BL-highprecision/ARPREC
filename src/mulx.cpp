/*
 * src/mprealx.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * This file contains the basic routines for very high levels of precision,
 * > 1000 digits.
 *
 * The main routines among these  mpmulx, which employs FFT for the
 * convolution in multiplication, and mpinix, which initializes certain
 * root of unity arrays for the fft.
 * All other routines in this file are here only to help mpmulx.
 *
 */
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpmulx(const mp_real& a, const mp_real& b, mp_real& c, int prec_words)
{
  /**
   * This routine multiplies MP number A and B to yield the MP product C.
   * Before calling MPMULX, the arrays UU1 and UU2 must be initialized by
   * calling MPINIX.  For modest levels of precision, use MPMUL.  Debug
   * output starts with debug_level == 8.
   *
   * This routine returns up to prec_words mantissa words of the product.  
   * If the complete double-long product of A and B is desired (for example
   *  in large integer applications), then prec_words must be at least as large
   * as the sum of the mantissa lengths of a and b.  In other words, if the
   * precision levels of a and b are both 256 words, then prec_words must be at least
   * 512 words to obtain the complete double-long product in C.
   *
   * This subroutine uses an advanced technique involving the FFT.
   * For high precision it is significantly faster than the conventional
   * scheme used in MPMUL.
   */
  

  double t1, t2;
  int ia, ib, na, nb, ncr, i, nn, nc, nx;
  int na2, nb2;

  if(error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if (debug_level >= 8) {
    cerr << "MPMULX I" << endl;
  }
  
  ia = sign(1, int(a[1]));
  ib = sign(1, int(b[1]));
  na = std::min(int(std::abs(a[1])), prec_words);
  nb = std::min(int(std::abs(b[1])), prec_words);
  ncr = 1 << mpmcrx; // pow(2.0, mpmcrx);
  if(!na || !nb) {
    // one of the inputs is zero -- result is zero.
    zero(c);
    return;
  }

  // Check if precision level of one of the arguments is 
  // enough to justify the advanced routine.
  
  if(na <= ncr || nb <= ncr) {
    mpmul(a, b, c, prec_words);
    return;
  }

  double *d1 = new double[4*prec_words+8];
  double *d2 = new double[4*prec_words+8];
  double *d3 = new double[8*prec_words+16];
  const double mp_down14 = 1.0 / 4096.0;
  const double mp_down24 = mprbx;
  const double mp_down34 = mp_down14 * mp_down24;
  const double mp_up14 = 4096.0;
  const double mp_up24 = mpbbx;
  const double mp_up34 = mp_up24 * mp_up14;

  // Place the input data in A and B into the scratch arrays DD1 and 
  // DD2.  this code alsp splits the input data into fourth sized words.
  // (mpnbt/4) bits for each new word (maximum)
  int i2=1;
  na2 = 4*na;
  nb2 = 4*nb;
  for(i=0;i<na;i++) {
    i2 = 4*i;
    t1 = a[i+FST_M];
    t2 = FLOOR_POSITIVE(mp_down34*t1);
    d1[i2] = t2;
    t1 -= mp_up34 * t2;
    t2 = int(mp_down24*t1);
    d1[i2+1] = t2;
    t1 -= mp_up24 * t2;
    t2 = int(mp_down14*t1);
    d1[i2+2] = t2;
    t1 -= mp_up14 * t2;
    d1[i2+3] = t1;
  }

  for(i = na2;i<nb2;i++)
    d1[i] = 0.0;

  for(i=0;i<nb;i++) {
    i2 = 4*i;
    t1 = b[i+FST_M];
    t2 = FLOOR_POSITIVE(mp_down34*t1);
    d2[i2] = t2;
    t1 -= mp_up34 * t2;
    t2 = int(mp_down24*t1);
    d2[i2+1] = t2;
    t1 -= mp_up24 * t2;
    t2 = int(mp_down14*t1);
    d2[i2+2] = t2;
    t1 -= mp_up14 * t2;
    d2[i2+3] = t1;
  }

  for(i = nb2;i<na2;i++)
    d2[i] = 0.0;

  // Call the convolution

  nn = std::max(na2, nb2);
  nx = int(ANINT(sqrt(3.0 * nn) + mprxx));
  mplconv(2, nn, nx, d1, d2, d3);

  // Recombine words and release carries.

  nc = std::min(na + nb, prec_words+3);
  int nc1 = nc-1;
  d1[1] = ia+ib ? nc : -nc;
  d1[2] = a[2] + b[2] + 1;

  d1[3] = d3[0] * mp_up24 + d3[1] * mp_up14 + d3[2];
  d1[nc+FST_M+1] = d1[nc+FST_M+2] = 0.0;
  
  for(i=0;i<nc1;i++) {
    i2 = i * 4 + 3;
    t1 = d3[i2];
    t2 = FLOOR_POSITIVE(t1 * mp_down14);
    t1 -= t2 * mp_up14;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] = t1 * mp_up34;
    t1 = d3[i2+1];
    t2 = int(t1 * mp_down24);
    t1 -= t2 * mp_up24;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] += t1 * mp_up24;
    t1 = d3[i2+2];
    t2 = int(t1 * mp_down34);
    t1 -= t2 * mp_up34;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] += t1 * mp_up14 + d3[i2+3];
  }
  //one more to do.
  t1 = d3[i2+4];
  t2 = FLOOR_POSITIVE(t1 * mp_down14);
  t1 -= t2 * mp_up14;
  d1[i+FST_M] += t2;
  d1[i+FST_M+1] = t1 * mp_up34;

  //now go on.
  int d_add=0;
  i=0;
  //eliminate leading zeros (except where something
  // will probably carry in)
  while(i<(nc1-3) && d1[i+FST_M] == 0.0 && d1[i+FST_M+1] < (mpbdx-2.0)) {
    i++;
  }
  if(i) {
    d1[2] -= double(i);
    d1[1] = sign(1.0, d1[1]) * (std::abs(d1[1])  - double(i));
    d1[i+2] = d1[2];
    d1[i+1] = d1[1];
    d1 += i;
    d_add -= i;
  }

  // Fix up the result;
  mpnorm(d1, c, prec_words);
  
  if (debug_level >= 8) cerr << "MPMULX 0" << endl;
  delete [] (d1+d_add);
  delete [] d2;
  delete [] d3;
  return;
}

