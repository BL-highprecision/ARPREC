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

void mp_real::mpnpwr(const mp_real& a, int n, mp_real& b, int prec_words)
{
  /**
   * This computes the N-th power of the MP number A and returns the MP result
   * in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
   * of A ^ |N| is returned.  For extra high levels of precision, use MPNPWX.
   * Debug output starts with debug_level = 7.
   *
   * This routine employs the binary method for exponentiation.
   */

  int kk, kn, n5, na, nn, nws, skip;
  //const double cl2 = 1.4426950408889633;  // log2(e)
  mp_real f1(1.0, 9), sk0(0.0, prec_words+6), sk1(0.0, prec_words+6);
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(b);
    return;
  }
  if (debug_level >= 7) print_mpreal("MPNPWR I ", a);
  
  na = std::min (int(std::abs(a[1])), prec_words);
  if (na == 0) {
    if (n >= 0) {
      zero(b);
      if (debug_level >= 7) print_mpreal("MPNPWR O ", b);
      return;
    } else {
      if (MPKER[57] != 0) {
        cerr << "*** MPNPWR: Argument is zero and N is negative or zero." << endl;
        error_no = 57;
        if (MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }

  n5 = prec_words + 5;
  nws = prec_words;
  prec_words = prec_words + 1;
  nn = std::abs (n);
  
  skip = 0;
  if (nn == 0) {
    mpeq(f1, b, prec_words);
    prec_words = nws;
    if (debug_level >= 7) print_mpreal("MPNPWR O ", b);
    return;
  } else if (nn == 1) {
    mpeq(a, b, prec_words);
    skip = 1;
  } else if (nn == 2) {
    mpmul(a, a, sk0, prec_words);
    mpeq(sk0, b, prec_words);
    skip = 1;
  }

  if ( skip == 0 ) {
    // nn has the power at beginning
    mpeq (a, sk0, prec_words);
    mpeq (f1, b, prec_words);
    kn = nn;
    
    // Compute B ^ N using the binary rule for exponentiation.
    while(kn) {
      kk = kn / 2;
      if (kn != 2 * kk) { // kn is odd
        mpmul(b, sk0, sk1, prec_words);
        mpeq(sk1, b, prec_words);
      }
      kn = kk;
      if (kn) {
        mpmul(sk0, sk0, sk1, prec_words);
        mpeq(sk1, sk0, prec_words);
      }
    }
  }
  
  // Compute reciprocal if N is negative.
  if (n < 0) {
    mpdiv(f1, b, sk0, prec_words);
    mpeq(sk0, b, prec_words);
  }
  
  // Restore original precision level.
  prec_words = nws;
  mproun(b);
  
  if (debug_level >= 7) print_mpreal("MPNPWR O ", b);

}

