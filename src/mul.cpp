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
#include <cassert>
#include <arprec/mp_real.h>
#include "small_inline.h"

void mp_real::mpmul(const mp_real& a, const mp_real&b, mp_real& c, int prec_words)
{
  /**
   * This routine multiplies MP numbers A and B to yield the MP product C.
   * When one of the arguments has a much higher level of precision than the
   * other, this routine is slightly more efficient if A has the lower level of
   * precision.  For extra high levels of precision, use MPMULX.  Debug output
   * starts with debug_level = 8.
   *
   * This routine returns up to MPNW mantissa words of the product.  If the
   * complete double-long product of A and B is desired (for example in large
   * integer applications), then MPNW must be at least as large as the sum of
   * the mantissa lengths of A and B.  In other words, if the precision levels
   * of A and B are both 64 words, then MPNW must be at least 128 words to
   * obtain the complete double-long product in C.
   */
  int i, j, j3, jd, ia, ib, na, nb, nc, n2;
  double d2, t1, t2, t[2], a_val;
  double* d;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(c);
    return;
  }
  if (debug_level >= 8) {
    print_mpreal("MPMUL a ", a);
    print_mpreal("MPMUL b ", b);
  }
  
  ia = int(sign(1.0, a[1]));
  ib = int(sign(1.0, b[1]));
  na = std::min (int(std::abs(a[1])), prec_words);
  nb = std::min (int(std::abs(b[1])), prec_words);
  
  // One of the inputs is zero -- result is zero.
  if (na == 0 || nb == 0) {
    zero(c);
    if (debug_level >= 8) print_mpreal("MPMUL O ", c);
    return;
  }

  // One of the inputs is 1 or -1.
  if(na == 1) {
    if (a[3] == 1.) {
      // A is 1 or -1 -- result is B or -B.
      nc = std::min(nb, (int(c[0])-FST_M-2));
      c[1] = sign(nc, ia * ib);
      c[2] = a[2] + b[2];
      
      for (i = FST_M; i < nc + FST_M; ++i) c[i] = b[i];
      if(nc < nb) {
        //round..
        c[FST_M+nc] = b[FST_M+nc];
        c[FST_M+nc+1] = 0.0;
        mproun(c);
      } else {
        c[FST_M+nc] = c[FST_M+nc+1] = 0.0;
      }    
      if (debug_level >= 8) print_mpreal("MPMUL O ", c);
      return;
    } else {
      //a has only one mantissa word, just use mpmuld
      double a2 = a[2];
      mpmuld(b, ia * a[FST_M], 0, c, prec_words);
      c[2] += a2;
      return;
    }
  } else if (nb == 1) {
    if(b[3] == 1.) {
      // B is 1 or -1 -- result is A or -A.
      nc = std::min(na, (int(c[0])-FST_M-2));
      c[1] = sign(nc, ia * ib);
      c[2] = a[2] + b[2];
      
      for (i = FST_M; i < nc + FST_M; ++i) c[i] = a[i];
      if(nc < na) {
        //round..
        c[FST_M+nc] = a[FST_M+nc];
        c[FST_M+nc+1] = 0.0;
        mproun(c);
      } else {
        c[FST_M+nc] = c[FST_M+nc+1] = 0.0;
      }
      
      if (debug_level >= 8) print_mpreal("MPMUL O ", c);
      return;
    } else {
      //b has only one mantissa word, just use mpmuld
      double b2 = b[2];
      mpmuld(a, ib * b[FST_M], 0, c, prec_words);
      c[2] += b2;
      return;
    }
  }
  // ok. ready for main part of routine
  d = new double[prec_words+5+FST_M];  // accumulator

  nc = std::min(int(c.mpr[0])-5, std::min (na + nb, prec_words));
  d2 = a[2] + b[2]; // exponent
  for (i = 1; i < prec_words+4+FST_M; ++i) d[i] = 0.0;


  // Perform ordinary long multiplication algorithm. Accumulate at most 
  //  MPNW+5-FST+1 == (max j possible)-FST_M+1
  // mantissa words of the product.

  for (j = FST_M; j < na + FST_M; ++j) {
    a_val = a[j];
    j3 = j - FST_M;
    n2 = std::min (nb + FST_M, prec_words + 5 - j3);
    
    jd = j;
    for(i = FST_M; i < n2; ++i) {
      t[0] = mp_two_prod_positive(a_val, b[i], t[1]); 
                // t[0], t[1] non-overlap, <= 2^mpnbt-1
      d[jd-1] += t[0];
      d[jd] += t[1];
      ++jd;
    }

      // Release carry to avoid overflowing the exact integer capacity
      // (2^mpnbt-1) of a floating point word in D.
    if(!((j-2) & (mp::mpnpr-1))) { // assume mpnpr is power of two
      for(i= jd-1;i>=j;i--) {
          t1 = d[i];
            t2 = int (t1 * mprdx);     // carry <= 1
            d[i] = t1 - t2 * mpbdx;   // remainder of t1 * 2^(-mpnbt)
          d[i-1] += t2;
      }
    }
  }
  int d_add = 0;
 
  // If D[1] is nonzero, shift the result two words right.
  if (d[1] != 0.0) {
    // this case shouldn't really happen.
    assert(0);
    d2 += 2.0;
    for (i = nc + 4; i >= FST_M; --i)
      d[i] = d[i-2];    
  } else if (d[2] != 0.0 || (d[3] >= mpbdx && (d[2] = 0.0, 1))) {
  // If D[2] is nonzero, shift the result one word right.
    d2 += 1.0;  // exponent
    d--; d_add++;
  }
  // Result is negative if one of {ia, ib} is negative.
  d[1] = ia+ib ? nc : -nc;
  d[2] = d2;

  //  Fix up result, since some words may be negative or exceed MPBDX.
  mpnorm(d, c, prec_words);
  delete [] (d +  d_add);
  
  if (debug_level >= 8) print_mpreal("MPMUL O ", c);
  return;
}

