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

void mp_real::mpmulacc(const mp_real &a, const mp_real &b, mp_real &c) {
  /**
   * This routine computes
   *
   *      C  <--  C  +  A * B
   *
   * If one of the multiplicand has a much higher level or precision than
   * the other, this routine is more efficient if A has the lower level of
   * precision.  
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
  
  // One of the inputs is zero -- the result is unchanged.
  if (na == 0 || nb == 0) {
    if (debug_level >= 8) print_mpreal("MPMUL O ", c);
    return;
  }

  // One of the inputs is 1 or -1.
  if(na == 1) {
    if (a[3] == 1.) {
      int ic = ia * ib;
      if (ic == 1)
        c += b;
      else
        c -= b;
      return;
    }
  } else if (nb == 1) {
    if(b[3] == 1.) {
      int ic = ia * ib;
      if (ic == 1)
        c += a;
      else
        c -= a;
      return;
    }
  }

  // ok. ready for main part of routine
  d = new double[prec_words+5+FST_M];  // accumulator

  nc = std::min(int(c.mpr[0])-5, std::min (na + nb, prec_words));
  d2 = a[2] + b[2]; // exponent
  //for (i = 1; i < prec_words+4+FST_M; ++i) d[i] = 0.0;

  int ic = (c[1] >= 0.0) ? 1 : -1;
  int iab = ia * ib;
  double c_exp = c[2];
  double ab_exp = d2;
  int delta;
  int n_c = std::abs(static_cast<int>(c[1]));
  /*
  cout << "ab_exp = " << ab_exp << endl;
  cout << "c_exp = " << c_exp << endl;
  cout << "n_c = " << n_c << endl;
  cout << "nc = " << nc << endl;
  cout << "na = " << na << endl;
  cout << "nb = " << nb << endl;
  cout << "prec_words = " << prec_words << endl;
  */
  d[0] = 0.0;
  if (ab_exp >= c_exp) {
    d[2] = d[1] = 0.0;
    delta = static_cast<int>(ab_exp - c_exp);
    for (i = 3; i < 3 + delta; i++) d[i] = 0.0;
    int lim = 3 + delta + n_c;
    lim = std::min(lim, prec_words+4+FST_M);
    if (ic == iab)
      for (i = 3 + delta; i < lim; i++) d[i] = c[i-delta];
    else
      for (i = 3 + delta; i < lim; i++) d[i] = -c[i-delta];
    for (i = lim; i < prec_words+4+FST_M; i++) d[i] = 0.0;
    /*
    for (i = 1; i < prec_words+4+FST_M; i++) 
      cout << "d[" << i << "] = " << d[i] << endl;
    for (i = 3; i < 3 + n_c; i++)
      cout << "c[" << i << "] = " << c[i] << endl;
      */
    delta = 0;
  } else {
    d[2] = d[1] = 0.0;
    delta = static_cast<int>(c_exp - ab_exp);
    int lim = std::min(n_c + 3, prec_words+4+FST_M);
    if (ic == iab)
      for (i = 3; i < lim; i++) d[i] = c[i];
    else
      for (i = 3; i < lim; i++) d[i] = -c[i];
    for (i = lim; i < prec_words+4+FST_M; i++) d[i] = 0.0;
    /*
    for (i = 1; i < prec_words+4+FST_M; i++) 
      cout << "d[" << i << "] = " << d[i] << endl;
    for (i = 3; i < 3 + n_c; i++)
      cout << "c[" << i << "] = " << c[i] << endl;
      */
    //c += a * b;
    //return;
  }
  //for (i = 1; i < prec_words+4+FST_M; ++i) d[i] = 0.0;

  // Perform ordinary long multiplication algorithm. Accumulate at most 
  //  MPNW+5-FST+1 == (max j possible)-FST_M+1
  // mantissa words of the product.

  int last = 1;
  for (j = FST_M; j < na + FST_M; ++j, ++last) {
    a_val = a[j];
    j3 = j - FST_M;
    n2 = std::min (nb + FST_M, prec_words + 5 - j3 - delta);
    
    jd = j + delta;
    for(i = FST_M; i < n2; ++i) {
      t[0] = mp_two_prod_positive(a_val, b[i], t[1]); 
                // t[0], t[1] non-overlap, <= 2^mpnbt-1
      d[jd-1] += t[0];
      d[jd] += t[1];
      ++jd;
    }

      // Release carry to avoid overflowing the exact integer capacity
      // (2^mpnbt-1) of a floating point word in D.
    if (last >= mp::mpnpr-2) {
      last = 1;
      for(i= jd-1;i>=j;i--) {
          t1 = d[i];
            t2 = int (t1 * mprdx);     // carry <= 1
            d[i] = t1 - t2 * mpbdx;   // remainder of t1 * 2^(-mpnbt)
          d[i-1] += t2;
      }
    }
  }
  int d_add = 0;
  /*
  for (i = 0; i < 3 + prec_words+4+FST_M; i++)
    cout << "d[" << i << "] = " << d[i] << endl;
    */

  //for (i = 1; i < prec_words+4+FST_M; ++i) d[i] = 0.0;
 
  // If D[1] is nonzero, shift the result two words right.
  //if (delta == 0) {
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
  //} else {
    d2 += delta;
  //}
  d[1] = (iab == -1) ? -nc : nc;
  //d[1] = ia+ib ? nc : -nc;
  // Result is negative if one of {ia, ib} is negative.
  d[2] = d2;

  /*
  for (i = 0; i < prec_words+4+FST_M; i++)
    cout << d[i] << " ";
  cout << endl;
  */
  //  Fix up result, since some words may be negative or exceed MPBDX.
  mpnorm(d, c, prec_words);
  /*
  for (i = 0; i < c[0]; i++)
    cout << c[i] << " ";
  cout << endl;
  */
  delete [] (d +  d_add);
  
  if (debug_level >= 8) print_mpreal("MPMUL O ", c);
  return;
}

