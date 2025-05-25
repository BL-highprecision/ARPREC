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

using std::cerr;
using std::endl;

void mp_real::mpinfr (const mp_real& a, mp_real& b, mp_real& c, 
                      int prec_words, int want_frac /*=1*/) 
{
  /*
   *  mpinfr sets B to the integer part of A, and C to the
   *  fractional part of A.  Note that if A = -3.3, then B = -3, 
   *  and C = -0.3.
   *  Debug output starts with debug_level == 9.
   *
   *  The input argument want_frac is optional.  It defaults to 1. if
   *  frational part is not desired, pass in 0 for want_frac.
   *
   *  Required space for B and C: prec_words+4 cells each.
   */
  double ma;
  int na, nb, nc=0, ia;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b); zero(c);
    return;
  }
  if(debug_level >= 9) cerr << "MPINFR I : A = " << a << endl;
  
  // Check if A is zero.
  ia = sign(1, int(a[1]));
  na = std::min(int(std::abs(a[1])), prec_words);
  ma = a[2];
  if(na == 0) {
    zero(b);
    if ( want_frac ) zero(c);
    return;
  }
  
  if(ma > prec_words - 1) {
    //The input is alreay an integer.  
    //Any fractional part has been lost due to insufficient precision.
    if(MPKER[40] != 0) {
      cerr << "***MPINFR: Argument is too large." << endl;
      error_no = 40;
      if(MPKER[error_no] == 2) mpabrt();
    }
    //don't return here, Since the user might want the answer anyway,
    //such as in mpoutx_helper.
  }
  
  // Place integer part in B.
  nb = std::min(std::max(int(ma)+1, 0), na);
  if(nb == 0) {
    zero(b);
  } else {
    int nb2 = std::min(nb, int(b[0])-FST_M-2);
    b[1] = sign(nb2, ia);
    b[2] = ma;
    b[nb2+FST_M] = b[nb2+FST_M+1] = 0.0;
    int nb3 = std::min(nb2+1, nb);
    for(int i = FST_M;i <= nb3 + 2; i++)
      b[i] = a[i];
  }
  
  // Place fractional part in C.
  
  if(want_frac) {
    nc = na - nb;
    if(nc <= 0) {
      zero(c);
    } else {
      int nc2 = std::min(nc, int(c[0])-FST_M-2);
      c[1] = sign(nc2, ia);
      c[2] = ma - nb;
      c[nc2+FST_M] = c[nc2+FST_M+1] = 0.0;
      int nc3 = std::min(nc2+1, nc);
      for(int i = FST_M; i <=  nc3+2; i++) {
        c[i] = a[i+nb];
      }
    }
  
    // Fix up results. B may have trailing zeros and C may have
    // Leading zeros.
    mproun(c);
  }
  // Fix up results. B may have trailing zeros and C may have
  // Leading zeros.
  mproun(b);
  
  if(debug_level >=9) cerr << "MPINFR 0 : nb = " << nb << ", nc = " << nc << endl;
} 
  
