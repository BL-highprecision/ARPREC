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

void mp_real::mpdiv(const mp_real& a, const mp_real& b, mp_real& c, int prec_words)
{
  /**
   * This divides the MP number A by the MP number B to yield the MP 
   * quotient C.  For extra high levels of precision, use MPDIVX.
   * Debug output starts with debug_level = 8.
   *
   * The algorithm is by long division.
   */

  int i, ia, ib, ij, is, i2, i3=0, j, j3, na, nb, nc, BreakLoop;
  double rb, ss, t0, t1, t2, t[2];
  double* d;
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(c);
    return;
  }
  
  if (debug_level >= 8) {
    print_mpreal("MPDIV a ", a);
    print_mpreal("MPDIV b ", b);
  }
  
  ia = (a[1] >= 0 ? 1 : -1); 
  ib = (b[1] >= 0 ? 1 : -1); 
  na = std::min (int(std::abs(a[1])), prec_words);
  nb = std::min (int(std::abs(b[1])), prec_words);
  
  //  Check if dividend is zero.
  if (na == 0) {
    zero(c);
    if (debug_level >= 8) print_mpreal("MPDIV O ", c);
    return;
  }
  
  if (nb == 1 && b[FST_M] == 1.) {
    // Divisor is 1 or -1 -- result is A or -A.
    c[1] = sign(na, ia * ib);
    c[2] = a[2] - b[2];
    for (i = FST_M; i < na+FST_M; ++i) c[i] = a[i];
    
    if (debug_level >= 8) print_mpreal("MPDIV O ", c);
    return;
  }
  
  //  Check if divisor is zero.
  if (nb == 0) {
    if (MPKER[31] != 0) {
      cerr << "*** MPDIV: Divisor is zero." << endl;
      error_no = 31;
      if (MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  //need the scratch space now...
  d = new double[prec_words+9];
  int d_add=0;
  d++; d_add--;

  // Initialize trial divisor and trial dividend.
  t0 = mpbdx * b[3];
  if (nb >= 2) t0 = t0 + b[4];
  if (nb >= 3) t0 = t0 + mprdx * b[5];
  rb = 1.0 / t0;
  d[0]  = d[1] = 0.0;
  
  for (i = 2; i < na+2; ++i) d[i] = a[i+1];
  for (/*i = na+2*/; i <= prec_words+7; ++i) d[i] = 0.0;

  // Perform ordinary long division algorithm.  First compute only the first
  // NA words of the quotient.
  for (j = 2; j <= na+1; ++j) {
    t1 = mpbx2 * d[j-1] + mpbdx * d[j] + d[j+1];
    t0 = AINT (rb * t1); // trial quotient, approx is ok.
    j3 = j - 3;
    i2 = std::min (nb, prec_words + 2 - j3) + 2;
    ij = i2 + j3;
    for (i = 3; i <= i2; ++i) {
      i3 = i + j3;
      t[0] = mp_two_prod(t0, b[i], t[1]);
      d[i3-1] -= t[0];   // >= -(2^mpnbt-1), <= 2^mpnbt-1
      d[i3] -= t[1];
    }
    
      // Release carry to avoid overflowing the exact integer capacity
      // (2^52-1) of a floating point word in D.
    if(!(j & (mp::mpnpr-1))) { // assume mpnpr is power of two
      t2 = 0.0;
      for(i=i3;i>j+1;i--) {
        t1 = t2 + d[i];
        t2 = int (t1 * mprdx);     // carry <= 1
        d[i] = t1 - t2 * mpbdx;   // remainder of t1 * 2^(-mpnbt)
      }
      d[i] += t2;
    }
    
    d[j] += mpbdx * d[j-1];
    d[j-1] = t0; // quotient
  }
  
  // Compute additional words of the quotient, as long as the remainder
  // is nonzero.  
  BreakLoop = 0;
  for (j = na+2; j <= prec_words+3; ++j) {
    t1 = mpbx2 * d[j-1] + mpbdx * d[j];
    if (j < prec_words + 3) t1 += d[j+1];
    t0 = AINT (rb * t1); // trial quotient, approx is ok.
    j3 = j - 3;
    i2 = std::min (nb, prec_words + 2 - j3) + 2;
    ij = i2 + j3;
    ss = 0.0;
    
    for (i = 3; i <= i2; ++i) {
      i3 = i + j3;
      t[0] = mp_two_prod(t0, b[i], t[1]);
      d[i3-1] -= t[0];   // >= -(2^mpnbt-1), <= 2^mpnbt-1
      d[i3] -= t[1];
      
      //square to avoid cancellation when d[i3] or d[i3-1] are negative
      ss += sqr (d[i3-1]) + sqr (d[i3]); 
    }
      // Release carry to avoid overflowing the exact integer capacity
      // (2^mpnbt-1) of a floating point word in D.
    if(!(j & (mp::mpnpr-1))) { // assume mpnpr is power of two
      t2 = 0.0;
      for(i=i3;i>j+1;i--) {
        t1 = t2 + d[i];
        t2 = int (t1 * mprdx);     // carry <= 1
        d[i] = t1 - t2 * mpbdx;   // remainder of t1 * 2^(-mpnbt)
      }
      d[i] += t2;
    }

    d[j] += mpbdx * d[j-1];
    d[j-1] = t0;
    if (ss == 0.0) {
      BreakLoop = 1;
      break;
    }
    if (ij <= prec_words+1) d[ij+3] = 0.0;
  } 

  // Set sign and exponent, and fix up result.
  if(!BreakLoop) j--;
  d[j] = 0.0;
  
  if (d[1] == 0.0) {
    is = 1;
    d--; d_add++;
  } else {
    is = 2;
    d-=2; d_add+=2;
    //for (i = j+1; i >= 3; --i) d[i] =  d[i-2];
  }

  nc = std::min( (int(c[0])-FST_M-2), std::min (j-1, prec_words));
  

  d[1] = ia+ib ? nc : -nc;//sign(nc, ia * ib);
  d[2] = a[2] - b[2] + is - 2;  
  
  mpnorm(d, c, prec_words);
  delete [] (d+d_add);
  
  if (debug_level >= 8) print_mpreal("MPDIV O ", c);
  return;
}

