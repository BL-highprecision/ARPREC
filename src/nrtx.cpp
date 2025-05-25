/*
 * src/mprealx3.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * This file contains additional high precision mp_real routines.
 */

#include <arprec/mp_real.h>

using std::cerr;
using std::endl;

void mp_real::mpnrtx(const mp_real& a, int n, mp_real& b)
{
  /**
   * This computes the N-th root of the MP number A and returns the MP result
   * in B.  N must be at least one and must not exceed 2 ^ 30.  For modest
   * levels of precision, use MPNRT.  Debug output starts with debug_level = 6.
   *
   * This routine uses basically the same Newton iteration algortihm as
   * MPNRT.  In fact, this routine calls mpnrt to obtain an initial approxi-
   * mation.  See the comment about the parameter NIT in MPDIVX.
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
   *
   *  It is important to note that the results of mpnrt are undefined
   * when a and b are the same object.
   *
   */
  int k, n1;
  double t2, tn;
  const double cl2 = 1.4426950408889633;
  const int nit = 3, n30 = 1 << 30; // 2 ^ 30
  int prec_words = mp::prec_words;
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  
  if (debug_level >= 6) print_mpreal("MPNRTX I: a ", a);
  
  int ncr = 1 << mpmcrx;
  int na = std::min (std::abs(static_cast<int>(a[1])), prec_words); // number of words in A
  
  if (na == 0) {
    zero(b);
    if (debug_level >= 6) print_mpreal("MPNRT O: b ", b);
    return;
  }
  
  if (a[1] < 0) {
    if (MPKER[59] != 0) {
      cerr << "*** MPNRT: Argument is negative.\n";
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
    mpeq (a, b, prec_words);
    if (debug_level >= 7) print_mpreal("MPNRTX O: b ", b);
    return;
  case 2:
    mpsqrtx (a, b, prec_words);
    if (debug_level >= 7) print_mpreal("MPNRTX O: b ", b);
    return;
#if 0
  case 3:
    mpcbrtx (a, b);
    if (debug_level >= 7) print_mpreal("MPNRTX O: b ", b);
    return;
#endif
  default: break;
  }

  // Check if precision level is too low to justify the advanced routine.

  if(prec_words  <= ncr) {
    mpnrt(a, n, b, prec_words); return;
  }

  int n5 = prec_words + 5;
  mp_real f1(1.0, 8), f2(0.0, 8), sk0(0.0, n5);
  mp_real sk1(0.0, n5), sk2(0.0, n5), sk3(0.0, n5);
  
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
        mpmulx(sk0, sk1, sk3, prec_words);
        mpeq(sk3, sk1, prec_words);
        mpadd(sk1, sk2, sk3, prec_words);
        mpeq(sk3, sk2, prec_words);
      } while (sk1[1] != 0 && sk1[2] >= -prec_words); 
      
      mpeq(sk2, b, prec_words);
      mpdivx(f1, sk2, sk0, prec_words);
      prec_words = nws;
      mproun(b);
      
      if (debug_level >= 7) print_mpreal("MPNRT O: b ", b);
      return;
    }
  }
  
  //  Compute the initial approximation of A ^ (-1/N).
  prec_words = ncr + 1;
  mpnrt(a, n, sk0, prec_words);
  mpdivx(f1, sk0, b, prec_words);

  tn = n;
  mpdmc (tn, 0, f2, prec_words);
  int iq = 0;
  
  //  Perform the Newton-Raphson iteration described above with a dynamically
  //  changing precision level MPNW (one greater than powers of two).
  for (k = mpmcrx+1; k <= mq; ++k) {
    prec_words = std::min (2 * prec_words - 2, nws) + 1;
    int loop = 1;
    while (loop) {
      mpnpwx (b, n, sk0, prec_words);
      mpmulx (a, sk0, sk1, prec_words);
      mpsub (f1, sk1, sk0, prec_words);
      mpmulx (b, sk0, sk1, prec_words);
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
  mpdivx (f1, b, sk1, prec_words);
  mpeq (sk1, b, prec_words);
  
  prec_words = nws;
  mproun (b);
  
  if (debug_level >= 6) print_mpreal("MPNRTX O: b ", b);
}

