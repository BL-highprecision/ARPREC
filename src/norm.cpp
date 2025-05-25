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

void mp_real::mpnorm(double d[], mp_real &a, int prec_words) {
  // This converts the MP number in array D of MPCOM4 to the standard
  // normalized form in A.  The MP routines often leave negative numbers or
  // values exceeding the radix MPBDX in result arrays, and this fixes them.
  // MPNORM assumes that two extra mantissa words are input at the end of D.
  // This reduces precision loss when it is necessary to shift the result to
  // the left.  This routine is not intended to be called directly by the user.
  // Debug output starts with debug_level = 10.
  //
  // Max SP space for A: MPNW + 5 cells.
  // 
  // The first 3 words of D have special meanings:
  //     d[0] : not accessed
  //     d[1] : sign and number of words in D
  //     d[2] : exponent
  // Before normalization, each word of D may have up to 53 bits and may be
  // negative. After normalization, each word is in [0, 2^50-1]
  //     
  double a2, t1, t2, t3;
  int i, ia, na, nd, n4;

  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(a);
    return;
  }

  if (debug_level >= 7) {
    mp_real t;
    t[1] = d[1]; t[2] = d[2];
    for (i = FST_M; i < std::abs(t[1]) + FST_M; ++i) t[i] = d[i];
    print_mpreal("MPNORM I ", t);
  }

  ia = static_cast<int>(sign(1.0, d[1]));
  nd = std::abs(static_cast<int>(d[1]));
  na = std::min(nd, prec_words);
  na = std::min(na, (static_cast<int>(a[0]) - 5)); // do not exceed the allocated memory
  if (na == 0) {
    zero(a);
    if (debug_level >= 9) print_mpreal("MPNORM O ", a);
    return;
  }
  n4 = na + 4; 
        /* n4 is often prec_words+FST_M+1. The carry release loop usually starts
         *  one word after the last word necessary for the specified
         * precision prec_words. */
        // if na == int(a[0]) - 5, then n4 is the last addressable word.
  a2 = d[2];

  // Carry release loop.
  t1 = 0.0;
  for (i = n4; i >= FST_M; --i) {
    t3 = t1 + d[i];  // exact 
    t2 = t3 * mprdx;
    t1 = int (t2);   // carry <= 1
    if ( t2 < 0.0 && t1 != t2 ) // make sure d[i] will be >= 0.0
      t1 -= 1.0;
    a[i] = t3 - t1 * mpbdx;
  }
  a[2] = t1;

  if ( a[2] < 0.0 ) {
    // The leading word is negative -- negate all words and re-normalize
    ia = -ia;
    
    for (i = 2; i <= n4; ++i)
      a[i] = -a[i];
    for (i = n4; i >= FST_M; --i) {
      if ( a[i] < 0 ) {
        a[i] = mpbdx + a[i]; // <= 2^50-1
        a[i-1] -= 1.0;
      }
    }
  }
  // Now either a[2]>0.0. or a[2] == 0.0.

  if ( a[2] > 0.0 ) {
    // A nonzero carry is "spilled" into a[2].  Shift the entire number
    // right one cell.  The exponent and length of the result are increased
    // by one.
    if (na != prec_words && na < static_cast<int>(a[0]) - 5) {
      // can afford to up by one word.
      for (i = n4+1; i >= FST_M; --i) a[i] = a[i-1];
      na = std::min (na+1, prec_words);
      a2 += 1;
    } else {
      //can't up by a word.  must truncate one word. 
      for (i = n4; i >= FST_M; --i) a[i] = a[i-1];
      a2 += 1;
    }
  }

  // Perform rounding and truncation.
  a[1] = ia >= 0 ? na : -na;
  a[2] = a2;

  if (debug_level >= 7) {
    print_mpreal("MPNORM before mproun ", a);
  }
  mproun(a);
  return;
}
