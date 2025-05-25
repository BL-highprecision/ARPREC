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

void mp_real::mpfft2(int is, int ns, int m, int n, double *x, double *y,
		     double* &out)
{
  /**
   * This performs NS simultaneous N-point complex-to-complex FFTs, 
   * where N = 2^M.  X is the input array, UU1 is the root
   * of unity array, and Y  is a scratch array.  X, Y, and UU1 are
   * double complex.  This routine is not intended to be called 
   * directly by the user. The output is put into either x or y,
   * depending on m, and the pointer OUT is set to equal x or y, 
   * whichever gets the output.
   *
   */
  
  int l;
  int m_odd=0;
  // dimensions : x[mpnrow+mpnsp1, n], y[mpnrow+mpnsp1, n].

  if(m & 0x1) {
    //m is odd, so we need to first perform one radix 2 iteration
    //before the radix 4 iterations start.

    mpfft3_radix2(is, m, ns, 0, 0, x, y);      
    for(l=3;l<=m;l+=4) {
      mpfft3(is, l, ns, m, n, y, x);
      if(l == m) {m_odd = 1; break;}
      mpfft3(is, l + 2, ns, m, n, x, y);
    }
    if(m_odd) {
      out = x;
    } else {
      out = y;
    }
    return;
  }
  //else , m^0x1 == 0

  for(l=2;l<=m;l+=4) {
    mpfft3(is, l, ns, 1, n, x, y);
    if(l == m) {m_odd = 1; break;}
    mpfft3(is, l + 2, ns, 1, n, y, x);
  }
  
  if(m_odd) {
    out = y;
  } else {
    out = x;
  }
  return;
}

