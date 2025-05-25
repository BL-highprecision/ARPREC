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

void mp_real::mpmuld(const mp_real& a, double b, int n, mp_real& c, int prec_words)
{
  /**
   * This routine multiplies the MP number A by the DPE number (B, N) to yield
   * the MP product C.  Debug output starts with debug_level = 9.
   *
   * Here, DPE means double precision + exponent, so (B, N) = b * 2^n
   */

  int i, ia, ib, k, na, n1, n2, d_add;
  double bb, t[2];
  double* d;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(c);
    return;
  }
  
  if (debug_level >= 8) {
    print_mpreal("MPMULD a ", a);
    cerr << " MPMULD b " << b << " " << n << endl; 
  }
  
  // Check for zero inputs.
  ia = int(sign(1.0, a[1]));
  na = std::min (int(std::abs(a[1])), prec_words);
  ib = int(sign(1.0, b));
  if (na == 0 || b == 0.0) {
    zero(c);
    if (debug_level >= 9) print_mpreal("MPMULD O ", c);
    return;
  }
  if(n) {
    n1 = n / mpnbt;      // n = mpnbt*n1+n2
    n2 = n - mpnbt * n1;
    bb = ldexp(fabs(b), n2);
  } else {
    n1 = n2 = 0;
    bb = fabs(b);
  }

  // Reduce BB to within 1 and MPBDX.
  if (bb >= mpbdx) {
    for (k = 1; k <= 100; ++k) {
      bb = mprdx * bb;
      if (bb < mpbdx) {
        n1 = n1 + k;
        break;
      }
    }
  } else if (bb < 1.0) {
    for (k = 1; k <= 100; ++k) {
      bb = mpbdx * bb;
      if (bb >= 1.0) {
        n1 = n1 - k;
        break;
      }
    }
  }
#if 0
  printf("mpmuld[1]: bb = %22.18e, na = %d\n", bb, na);
#endif
  // BB is now between 1 and MPBDX (and positive)
  // If BB cannot be represented exactly in a single mantissa word, use MPMUL.
  if (bb != floor(bb)) {
    mp_real f(0.0, 9);
    mpdmc(b, n, f, prec_words);
    mpmul(f, a, c, prec_words);
    if (debug_level >= 9) print_mpreal("MPMULD O ", c);
    return; 
  }
  
  d = new double[prec_words+6];
  d_add = 0;

  // Perform short multiply operation.
  d[2] = 0.;
  for (i = FST_M; i < na + FST_M; ++i) {
#if 0
    printf("mpmuld[2]: a[%d] = %22.18e, bb = %22.18e\n", i, a[i], bb);
#endif
    t[0] = mp_two_prod_positive(a[i], bb, t[1]); 
        // t[0], t[1] non-overlap, <= 2^mpnbt-1
#if 0
    printf("mpmuld[3]: t[0] = %22.18e, t[1] = %22.18e\n", t[0], t[1]);
#endif
    d[i-1] += t[0];  // exact, <= 2^53-2
    d[i] = t[1];
#if 0
    printf("mpmuld[4]: d[%d] = %22.18e, d[%d] = %22.18e\n",
	   i-1, d[i-1], i, d[i]);
#endif
  }
  
  // If carry is nonzero, shift the result one word right.
  if (d[2] != 0.0) {
    ++n1;  // exponent
    ++na;//number of words

    d--;d_add++;// "shift" the array one to the right.
    // This has the same effect as the following commented out loop:
    //for (i = na + FST_M; i >= FST_M; --i) d[i] = d[i-1];
  }

  // Set the exponent and fix up the result.
  d[1] = ia+ib ? na : -na;//same as sign (na, ia * ib);
  d[2] = a[2] + n1;
  d[na+3] = 0.0;
  d[na+4] = 0.0;
  
  //  Fix up result, since some words may be negative or exceed MPBDX.
  mpnorm(d, c, prec_words);
  delete [] (d + d_add);

  if (debug_level >= 8) print_mpreal("MPMULD O ", c);
  return;
}

