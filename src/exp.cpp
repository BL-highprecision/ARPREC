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

void mp_real::mpexp(const mp_real& a, const mp_real& al2, mp_real& b, int prec_words)
{
  /**
   * This computes the exponential function of the MP number a and
   * returns the MP result in B.  al2 is the MP value of log(2),
   * produced by a prior call to MPLOG.  for extra high levels of precision,
   * use MPEXPX.  The last word of the result is not reliable.
   * Debug output starts with debug_level = 7. 
   *  
   * This routine uses a modification of the Taylor's series for Exp(t):
   * 
   * Exp(t) = (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q  *  2^n
   * 
   * where q = 256, r = t' / q, t' = t - n Log(2), and where n is 
   * chosen so that -0.5 Log(2) < t' <= 0.5 Log(2).  Reducing t mod Log(2)
   * and dividing by 256 insures that -0.001 < r <=0.001, which accelerates
   * convergence in the above series.
   *
   * Note that as the taylor series progresses, each successive term
   * is smaller in magnitude. This means that few mantissa words
   * of precision are needed for the later taylor series terms.
   * the computation of the taylor's series takes this into account 
   * (see below - term_prec), reducing the total amount of computation
   * by nearly half.
   */
  
  const double alt = 0.693147180559945309;
  const int nq = 8;
  double t1, t2;
  int n1, n2;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 7) cerr << "MPEXP I "<< a << endl;
  
  mpmdc(a, t1, n1, prec_words);
  if(n1 > 25) { // argument has very large magnitude
    if(t1 > 0.0) {
      if(MPKER[35] != 0) {
        cerr <<"*** MPEXP : Argument is too large : "<< a << endl;
        error_no = 35;
        if(MPKER[error_no] == 2) mpabrt();
      }
    } else {
      // large negative, just return 0.0
      zero(b);
    }
    return;
  }
  if(n1 != 0)
    t1 = ldexp(t1, n1);
  
  // Unless the argument is near log(2), log(2) must be precomputed. This
  // Exception is necessary because MPLOG calls MPEXP to initialize Log(2).
  
  if(std::abs(t1 - alt) > mprdx) {
    mpmdc(al2, t2, n2, prec_words);
    if(n2 != -mpnbt || (std::abs(t2 * mprdx - alt) > mprx2)) {
      if(MPKER[34] != 0) {
        cerr << "*** MPEXP: LOG (2) must be precomputed." << endl;
        error_no = 34;
        if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }
  
  // Check for overflows and underflows.
  // The constant is roughly (2^26 * mpnbt)* log(2),
  // Which allows a result with (word) exponent up to 2^26.
  if(std::abs(t1) > 2325815993.0) {
    if(t1 > 0.0) {
      if(MPKER[35] != 0) {
        cerr << "*** MPEXP : Argument is too large : " << a << endl;
        error_no = 35;
        if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    } else {
      // argument is very negative, the result would be
      // too small, so just return zero.
      zero(b);
      return;
    }
  }
  
  int n6 = prec_words + 6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6), sk3(0.0, n6);
  int nws = prec_words, nz, tl;
  prec_words++;
  mp_real f(1.0, 8);
  
  //        Compute the reduced argument A' = A - Log(2) * nint(A / Log(2)).
  //        Save NZ = nint(A / Log(2)) for correcting the exponent of the final
  //         result.

  if(std::abs(t1 - alt) > mprdx) {
    //It is unnecessary to compute much of fractional part of the following
    //division.
    prec_words = std::min(prec_words, int(a[2]) - int(al2[2]) + 3);
    prec_words = std::max(prec_words, 1);
    mpdiv(a, al2, sk0, prec_words);
    prec_words = nws+1;
    mpnint(sk0, sk1, prec_words);
    mpmdc(sk1, t1, n1, prec_words);

    nz = int(ldexp(t1, n1) + sign(mprxx, t1));
    mpmul(sk1, al2, sk2, prec_words);
    mpsub(a, sk2, sk0, prec_words);
  } else {
    mpeq(a, sk0, prec_words);
    nz = 0;
  }
  
  tl = int(sk0[2]) - prec_words;


  // Check if the reduced argument is zero

  if(sk0[1] == 0.0) {
    sk0[1] = 1.0;
    sk0[2] = 0.0;
    sk0[3] = 1.0;
    mpmuld(sk0, 1.0, nz, sk1, prec_words);
    mpeq(sk1, b, prec_words);
    prec_words = nws;
    mproun(b);
    return;
  }

  // Divide the reduced argument by 2^nq.

  mpdivd(sk0, 1.0, nq, sk1, prec_words);
  
  // Compute Exp using the usual Taylor series.
  
  mpeq(f, sk2, prec_words);
  mpeq(f, sk3, prec_words);
  const int max_iteration_count = 10000;
  int l1 = 0; //iteration number.
  int not_there_yet = 1;
  int term_prec;
  int i;

  while(not_there_yet && l1 < max_iteration_count) {
    l1++;
    t2 = l1;
    //get term precision from exponents of components of term, subtracting
    // exponent of current sum
    term_prec = std::min(nws+1, nws+int(sk2[2]+sk1[2]-sk3[2])+2);
    term_prec = std::max(term_prec, 0);
    if(term_prec <= 0) {
      prec_words = nws+1;
      break;
    }
    prec_words = term_prec;
    mpmul(sk2, sk1, sk0, prec_words);
    mpdivd(sk0, t2, 0, sk2, prec_words);
    prec_words = nws+1; // full precision to add term in.
    mpadd(sk3, sk2, sk3, prec_words);
    //the above line needs relies on mpadd being safe 
    //for use when the first and third arguments are the same object.

    // Check for convergence of the series.
    if(sk2[1] != 0.0 && sk2[2] > tl) {
      //keep going.
    } else {
      not_there_yet = 0;
    }
  }
  //check if exceeded iteration bound.
  if(l1 > max_iteration_count) {
    if(MPKER[36] != 0) {
      cerr <<"*** MPEXP: Iteration limit exceeded." << endl;
      error_no = 36;
      if(MPKER[error_no] == 2) mpabrt();
      prec_words = nws;
      return;
    }
  }
  //Result of taylor series stored in sk3.

  //         Raise to the (2^NQ)-th power.

  for(i=0;i<nq;i++) {
    mpmul(sk3, sk3, sk3, prec_words);
  }
  
  // Multiply by 2^NZ.
  if(nz) {
    mpmuld(sk3, 1.0, nz, b, prec_words);
  } else {
    mpeq(sk3, b, prec_words);
  }
  //restore original precision level.
  prec_words = nws;
  mproun(b);
  
  return;
}
