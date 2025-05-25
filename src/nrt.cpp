/*
 * src/mpreal2.cc
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

void mp_real::mpnrt(const mp_real& a, int n, mp_real& b, int prec_words)
{
  /**
   * This computes the N-th root of the MP number A and returns the MP result
   * in B.  N must be at least one and must not exceed 2 ^ 30.  For extra high
   * levels of precision, use MPNRTX.  Debug output starts with debug_level = 7.
   *
   * Max SP space for B: MPNW + 5 cells.
   *
   * This subroutine employs the following Newton-Raphson iteration, which
   * converges to A ^ (-1/N):
   *
   * X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)
   *
   * The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
   * These iterations are performed with a maximum precision level MPNW that
   * is dynamically changed, approximately doubling with each iteration.
   * See the comment about the parameter NIT in MPDIVX.
   *
   * When N is large and A is very near one, the following binomial series is
   * employed instead of the Newton scheme:
   *
   * (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...
   * See the comment about the parameter NIT in MPDIVX.
   *
   *  It is important to note that the results of mpnrt are undefined
   * when a and b are the same object.
   */
  int k, n1;
  double t2, tn;
  const double alt = 0.693147180559945309, cl2 = 1.4426950408889633;
  const int nit = 3, n30 = 1 << 30; // 2 ^ 30
  int n5 = prec_words + 5;
  mp_real f1(1.0, 8), f2(2.0, 8), sk0(0.0, n5);
  mp_real sk1(0.0, n5), sk2(0.0, n5), sk3(0.0, n5);
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  
  if (debug_level >= 7) print_mpreal("MPNRT I: a ", a);
  
  int na = std::min(std::abs(static_cast<int>(a[1])), prec_words); // number of words in A
  
  if (na == 0) {
    zero(b);
    if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
    return;
  }
  
  if (a[1] < 0) {
    if (MPKER[59] != 0) {
      cerr << "*** MPNRT: Argument is negative." << endl;
      error_no = 59;
      if (MPKER[error_no] == 2)  mpabrt();
    }
    return;
  }
  
  if (n <= 0 || n > n30) {
    if (MPKER[60] != 0) {
      cerr << "*** MPNRT: Improper value of N: " << n << endl;
      error_no = 60;
      if (MPKER[error_no] == 2)  mpabrt();
    }
    return;
  }
  
  //  If N = 1, 2 or 3,  MPEQ, MPSQRT or MPCBRT.  These are faster.
  switch (n) {
  case 1:
    mpeq(a, b, prec_words);
    if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
    return;
  case 2:
    mpsqrt (a, b, prec_words);
    if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
    return;
#if 0
  case 3:
    mpcbrt (a, b);
    if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
    return;
#endif
  default: break;
  }
  
  int nws = prec_words;
  
  //  Determine the least integer MQ such that 2 ^ MQ >= MPNW.
  double t1;
  t1 = prec_words;
  int mq = int(cl2 * log(t1) + 1.0 - mprxx); // *cast*
  
  //  Check how close A is to 1.
  mpsub(a, f1, sk0, prec_words);
  if (sk0[1] == 0) {
    mpeq(f1, b, prec_words);
    if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
    return;
  }
  
  mpmdc(sk0, t1, n1, prec_words);
  int n2 = int(cl2 * log (fabs (t1))); // *cast*
  t1 = ldexp(t1, -n2);
  n1 = n1 + n2;
  if (n1 <= -30) {
    t2 = n;
    n2 = int(cl2 * log (t2) + 1.0 + mprxx); // *cast*
    int n3 = - mpnbt * prec_words / n1;
    if (n3 < 1.25 * n2) {
      //  A is so close to 1 that it is cheaper to use the binomial series.
      prec_words = prec_words + 1;
      mpdivd(sk0, t2, 0, sk1, prec_words);
      mpadd(f1, sk1, sk2, prec_words);
      k = 0;
      
      do {
        k = k + 1;
        t1 = 1 - k * n;
        t2 = (k + 1) * n;
        mpmuld(sk1, t1, 0, sk3, prec_words);
        mpdivd(sk3, t2, 0, sk1, prec_words);
        mpmul(sk0, sk1, sk3, prec_words);
        mpeq(sk3, sk1, prec_words);
        mpadd(sk1, sk2, sk3, prec_words);
        mpeq(sk3, sk2, prec_words);
      } while (sk1[1] != 0 && sk1[2] >= -prec_words); 
      
      mpeq(sk2, b, prec_words);
      mpdiv(f1, sk2, sk0, prec_words);
      prec_words = nws;
      if (sk2[1] > b[1])
        mproun(b);
      
      if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
      return;
    }
  }
  
  //  Compute the initial approximation of A ^ (-1/N).
  tn = n;
  mpmdc (a, t1, n1, prec_words);
  n2 = int(-n1 / tn); // *cast*
  t2 = exp (-1.0 / tn * (log (t1) + (n1 + tn * n2) * alt));
  mpdmc (t2, n2, b, prec_words);
  mpdmc (tn, 0, f2, prec_words);
  prec_words = 3;
  int iq = 0;
  
  //  Perform the Newton-Raphson iteration described above with a dynamically
  //  changing precision level MPNW (one greater than powers of two).
  for (k = 1; k <= mq; ++k) {
    prec_words = std::min (2 * prec_words - 2, nws) + 1;
    int loop = 1;
    while (loop) {
      mpnpwr (b, n, sk0, prec_words);
      mpmul (a, sk0, sk1, prec_words);
      mpsub (f1, sk1, sk0, prec_words);
      mpmul (b, sk0, sk1, prec_words);
      mpdivd (sk1, tn, 0, sk0, prec_words);
      mpadd (b, sk0, sk1, prec_words);
      mpeq (sk1, b, prec_words);
      if (k == mq - nit && iq == 0) 
        iq = 1;
      else 
        loop = 0;
    }
  }
  //  Take the reciprocal to give final result.
  mpdiv (f1, b, sk1, prec_words);
  mpeq (sk1, b, prec_words);
  
  prec_words = nws;
  if (sk1[1] > b[1])
    mproun (b);
  
  if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
}

