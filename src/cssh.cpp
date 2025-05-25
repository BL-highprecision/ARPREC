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

void mp_real::mpcssh(const mp_real& a, 
                     const mp_real& al2, /* log of 2 */
                     mp_real& x,
                     mp_real& y)
{
  /**
   * This computes the hyperbolic cosine and sine of the
   * argument a, and returns the two MP results in x and y,
   * respectively.  AL2 is the MP value of Log (2) computed
   * by a previous call to MPLOG.  For extra high levels of precision,
   * use MPCSHX.  The last word of the result is not reliable.
   * Debug starts with debug_level == 5.
   * 
   *
   */
  int prec_words = mp::prec_words;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(x); zero(y); return;
  }
  
  int nws = prec_words;
  prec_words++;
  int n6 = prec_words+6; 
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6), sk3(0.0, n6), f(1.0, 8);

  mpexp(a, al2, sk0, prec_words);
  mpdiv(f, sk0, sk1, prec_words);
  mpadd(sk0, sk1, sk2, prec_words);
  mpmuld(sk2, 0.5, 0, sk3, prec_words);
  mpeq(sk3, x, prec_words);
  mpsub(sk0, sk1, sk2, prec_words);
  mpmuld(sk2, 0.5, 0, sk3, prec_words);
  mpeq(sk3, y, prec_words);
  
  //Restore original precision level.
  
  prec_words = nws;
  mproun(x);
  mproun(y);
}


