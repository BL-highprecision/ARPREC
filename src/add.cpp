/*
 * src/add.cc
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

void mp_real::mpadd(const mp_real &a, const mp_real &b, mp_real& c, int prec_words)
{
  /**
   * This routine adds MP numbers A and B to yield the MP sum C. It attempts
   * to include all significance of A and B in the result, up to the maximum
   * mantissa length MPNW.  Debug output starts with debug_level = 9. 
   * This is a new simplified version.
   *
   * D contains the intermediate results of addition before normalization.
   * The first 3 words of D have special meanings:
   *     d[0] : not accessed
   *     d[1] : sign and number of words in D
   *     d[2] : exponent
   * Before normalization, each word of D may have up to 53 bits and may be
   * negative. After normalization, each word is in [0, 2^mpnbt-1]
   */
  double db;
  int i, ia, ib, ish, ixa, ixb, ixd, na, nb;
  int m1, m2, m3, m4, m5, nsh;
  double *d;
  int nd; // number of actual words in d[]

  if (error_no != 0) {
    zero(c);
    return;
  }

  if (debug_level >= 9) {
    print_mpreal("MPADD a ", a);
    print_mpreal("MPADD b ", b);
  }
  
  ia = a[1] >= 0 ? 1 : -1;
  ib = b[1] >= 0 ? 1 : -1;
  na = std::min (static_cast<int>(std::abs(a[1])), prec_words); // number of words in A
  nb = std::min (static_cast<int>(std::abs(b[1])), prec_words); // number of words in B

  /* Check for zero inputs. */

  if (na == 0) {
    /* A is zero -- the result is B. */
    int num_words = std::min(nb, int(c[0])-FST_M);
    c[1] = ib > 0 ? num_words : -num_words;
    for (i = 2; i < num_words + FST_M; ++i) c[i] = b[i];
    return;
  } else if (nb == 0) {
    /* B is zero -- the result is A. */
    int num_words = std::min(na, int(c[0])-FST_M);
    c[1] = ia >= 0 ? num_words : -num_words; 
    for (i = 2; i < num_words + FST_M; ++i) c[i] = a[i];
    return;
  }

  // get ready for main part of routine.
  d = new double[prec_words+7];

  if (ia == ib) db = 1.0; //same signs - add
  else db = -1.0; // different signs - subtract

  ixa = static_cast<int>(a[2]);
  ixb = static_cast<int>(b[2]);
  ish = ixa - ixb;

  d[1] = 0.0;
  d[2] = 0.0;

  if (ish >= 0) { // |A| >= |B|
    // A has greater exponent than B, so B must be shifted to the right
    // to line up the radix point.
    
    m1 = std::min (na, ish);
    m2 = std::min (na, nb + ish);
    m3 = na;
    m4 = std::min (std::max (na, ish), prec_words + 1);
    m5 = std::min (std::max (na, nb + ish), prec_words + 1);
    //assert(m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5);
    
    for (i = FST_M; i < m1 + FST_M; ++i)
      d[i] = a[i];
    
    if(db > 0) {//Addition
      for (i = m1 + FST_M; i < m2 + FST_M; ++i)
        d[i] = a[i] + b[i-ish];
      
      for (i = m2 + FST_M; i < m3 + FST_M; ++i)
        d[i] = a[i];
    
      for (i = m3 + FST_M; i < m4 + FST_M; ++i)
        d[i] = 0.0;
      
      for (i = m4 + FST_M; i < m5 + FST_M; ++i)
        d[i] = b[i-ish];
    } else {//Subtraction
      for (i = m1 + FST_M; i < m2 + FST_M; ++i)
        d[i] = a[i] - b[i-ish];
      
      for (i = m2 + FST_M; i < m3 + FST_M; ++i)
        d[i] = a[i];
    
      for (i = m3 + FST_M; i < m4 + FST_M; ++i)
        d[i] = 0.0;
      
      for (i = m4 + FST_M; i < m5 + FST_M; ++i)
        d[i] = - b[i-ish];
    }
    nd = m5;
    ixd = ixa;
    d[nd+3] = 0.0;
    d[nd+4] = 0.0;

  } else {
    // B has greater exponent than A, so A must be shifted to the right
    // to line up the radix point.
    
    nsh = -ish;
    m1 = std::min (nb, nsh);
    m2 = std::min (nb, na + nsh);
    m3 = nb;
    m4 = std::min (std::max (nb, nsh), prec_words + 1);
    m5 = std::min (std::max (nb, na + nsh), prec_words + 1);
    //assert(m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5);
    
    if(db > 0) {//Addition
      for (i = FST_M; i < m1 + FST_M; ++i)
        d[i] = b[i];
      
      for (i = m1 + FST_M; i < m2 + FST_M; ++i)
        d[i] = a[i-nsh] + b[i];
      
      for (i = m2 + FST_M; i < m3 + FST_M; ++i)
        d[i] = b[i];

    } else {//Subtraction
      for (i = FST_M; i < m1 + FST_M; ++i)
        d[i] = - b[i];
      
      for (i = m1 + FST_M; i < m2 + FST_M; ++i)
        d[i] = a[i-nsh]  - b[i];
      
      for (i = m2 + FST_M; i < m3 + FST_M; ++i)
        d[i] = - b[i];
    }

    for (i = m3 + FST_M; i < m4 + FST_M; ++i)
      d[i] = 0.0;
    
    for (i = m4 + FST_M; i < m5 + FST_M; ++i)
      d[i] = a[i-nsh];
    
    nd = m5;
    ixd = ixb;
    d[nd+3] = 0.0;
    d[nd+4] = 0.0;
  }
  

  // Call mpnorm to fix up result and store in c.
  d[1] = ia >= 0 ? nd : -nd;
  d[2] = ixd;
  mpnorm(d, c, prec_words);

  delete [] (d);

  if (debug_level >= 9) print_mpreal("MPADD O: c ", c);
  return;
}

