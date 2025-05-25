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

void mp_real::mpdivd(const mp_real& a, double b, int n, mp_real& c, int prec_words)
{
  /**
   * This routine divides the MP number A by the DPE number (B, N) to yield
   * the MP quotient C.   Debug output starts with debug_level = 9.
   *
   */
  int ia, ib, i, j, k, na, nc, n1, n2, d_add;
  double bb, br, dd, t0, t[2];
  double* d;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(c);
    return;
  }
  if (debug_level >= 9) {
    print_mpreal("MPDIVD a ", a);
    cerr << " MIDIVD b " << b << endl;
    cerr << "        n " << n << endl;
  }

  ia = (a[1] >= 0 ? 1 : -1);
  na = std::min (int(std::abs(a[1])), prec_words);
  ib = (b >= 0 ? 1 : -1);
    
  // Check if divisor is zero.
  if (b == 0.0) {
    if (MPKER[32] != 0) {
      cerr << "*** MPDIVD: Divisor is zero.\b.n";     
      error_no = 32;
      if (MPKER[error_no] == 2) mpabrt();
    }
    return;
  }


  // Check if dividend is zero.
  if (na == 0) {
    zero(c);
    if (debug_level >= 9) print_mpreal("MPDIVD O ", c);
    return;
  }

  if(n) {
    n1 = n / mpnbt;
    n2 = n - mpnbt * n1;   // n = mpnbt*n1+n2
    bb = ldexp(fabs(b), n2);
  } else {
    n1 = n2 = 0;
    bb = fabs(b);
  }

  //  Reduce BB to within 1 and MPBDX.
  if (bb >= mpbdx) {
    for (k = 1; k <= 100; ++k) {
      bb = mprdx * bb;
      if (bb < mpbdx) {
        n1 += k;
        break;
      }
    }
  } else if (bb < 1.0) {
    for (k = 1; k <= 100; ++k) {
      bb = mpbdx * bb;
      if (bb >= 1.0) {
        n1 -= k;
        break;
      }
    }
  }
  
  // If B cannot be represented exactly in a single mantissa word, use MPDIV.
  if (bb != floor(bb)) {
    mp_real f(0.0, 9);
    bb = sign(bb, b);
    mpdmc(bb , n1*mpnbt, f, prec_words);
    mpdiv(a, f, c, prec_words);
    if (debug_level >= 9) print_mpreal("MPDIVD O ", c);
    return;
  }
  
  //Allocate scratch space.
  d  = new double[prec_words+6];  
  d_add = 0;

  br = 1.0 / bb;
  for (i = FST_M; i < na + FST_M; ++i) d[i] = a[i];
  for (/*i = na+FST_M*/; i <= prec_words+FST_M+2 ; i++) d[i] = 0.0;
  d[2] = 0.;

  // Perform short division (not vectorizable at present).
  // Continue as long as the remainder remains nonzero.
  for (j = FST_M; j <= prec_words+FST_M+1; ++j) {
    dd = mpbdx * d[j-1];
    if (j < na + FST_M)
      dd += d[j];
    else {
      if (dd == 0.0) {
               break;
      }
    }
    t0 = AINT (br * dd); // [0, 2^mpnbt-1], trial quotient.
    t[0] = mp_two_prod(t0, bb, t[1]); // t[0], t[1] non-overlap, <= 2^mpnbt-1
    d[j-1] -= t[0];
    d[j] -= t[1];
    d[j] += mpbdx * d[j-1];

    d[j-1] = t0; // quotient
  }
  
  //  Set sign and exponent of result.
  j--;
  if(AINT(d[j] * br)  != 0.0) {
    if(AINT(d[j] * br) >= 0.0)
      d[j-1] += 1.0;
    else
      d[j-1] -= 1.0;
  }

  d[j] = 0.;
  if ( d[2] != 0. ) {
    --n1;
    d--;d_add++;
    j++;
    //roughly equivelent slower version is commented out:
    //for (i = j; i >= FST_M; --i) d[i] = d[i-1];
  }
  
  nc = j-FST_M;
  //Quotient result is negative if exactly one of {ia, ib} is -1
  d[1] = ia+ib ? nc : -nc; 
  d[2] = a[2] - n1 -1;
  
  mpnorm(d, c, prec_words);
  delete [] (d+d_add);
  if (debug_level >= 9) print_mpreal("MPDIVD O ", c);
  return; 
}

