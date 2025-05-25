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

void mp_real::mpsqrt(const mp_real& a, mp_real& b, int prec_words)
{
  /**
   * This computes the square root of the MP number A and places
   * the result in B.  For extra high levels of precision,
   * use MPSQRTX. Debug output starts with debug_level =7.
   *  
   * Space suggested for B: prec_words+FST_M+2.
   *  
   * This subroutine employs the following Newton-Raphson
   * iteration, which converges to Sqrt(a):
   *
   *  X_{k+1} = X_k + 0.5 * (A - X_k^2) / X_k;
   *  
   * where the division () / X_k is performed with half normal
   * precision.  These iterations are performed with a maximum
   * precsion level of prec_words, with the precision dynamically changing,
   * doubling for each iteration.
   *
   * The precision doubling technique is efficient, but can cause
   * errors to build up in trailing mantissa words.  This error
   * can be controlled by repeating one of the iterations.
   * Iteration that is repieated is controlled by the input parameter
   * nit.  if nit == 0, the last iteration is repeated.  This
   * is the most effective, but the most expensive. if nit == 1, 
   * the second to last iteration is repeated instead, etc.
   *
   * Note that square roots of numbers with large
   * exponent magnitude, (base 2 exponent > 2^31), square roots
   * will be incorrect.  at present the library does not support
   * such large (or small) numbers.
   *
   * IMPORTANT - please note that mpsqrt is not safe for
   * use with the input variable the same as the output variable.
   * never call with mpsqrt(a, a).  This will lead to incorrect results.
   */
  const double cl2 = 1.4426950408889633; // == log_2(e) ==  1/log(2)
  const int nit = 3;

  if(error_no != 0) {
    if(error_no == 00) mpabrt();
    zero(b);
    return; 
  }
  if(debug_level >= 7 ) {
    cerr << "Runnung mpsqrt" << endl;
  }
  
  int ia = sign(1, int(a[1]));
  int na = std::min(int(std::abs(a[1])), prec_words);
  
  if(na == 0) {
    zero(b);
    return;
  }
  if(ia < 0.0) {//negative radicand!
    if(MPKER[70] != 0) {
      cerr << "*** MPSQRT: Argument is negative." << endl;
      error_no = 70;
      if(MPKER[error_no] == 2)
        mpabrt();
    }
    return;
  } // end negative radicand check

  int nws = prec_words;
  int k, mq, n, n2, iq=0;
  double t1, t2;
  int nw1, nw2, prec_change=0;
  int n7 = prec_words+7;
  mp_real sk0(0.0, n7), sk1(0.0, n7);

  // Determine the least integer MQ such that 2 ^ MQ >= prec_words.

  t1 = prec_words;
  mq = int(cl2 * log(t1) + 1.0 - mp::mprxx);
  
  //        Compute the initial approximation of Sqrt(A) using double
  //        precision.
#if 0
  printf("mpsqrt[1]  a= \n");
  cout << a;
  printf("prec_words %d\tt1 %22.18e\n", prec_words, t1);
#endif

  mpmdc(a, t1, n, prec_words);
#if 0
  printf("mpsqrt[2]  a= \n");
  cout << a;
  printf("prec_words %d\tt1 %22.18e\n", prec_words, t1);
#endif
  n2 = n / 2;
  t2 = sqrt((n2*2 == n) ? t1 : t1 * 2.0);
  t1 = t2;
  mpdmc(t1, n2, b, prec_words);
#if 0
  printf("mpsqrt[3] b=\n");
  cout << b;
  print_mpreal("comp cout ", b);
  printf("n2 %d\tt1 %22.18e\n", n2, t1);
#endif
  nw1 = nw2 = prec_words = 3;
  iq = 0;
  
  //         Perform the Newton-Raphson iteration described above by
  //         changing the precision level MPNW (one greater than powers of two).
  for(k=2;k <= mq;k++) {
    if(prec_change) {
      nw1 = prec_words;
      prec_words = std::min(2*prec_words-2, nws)+1; 
      nw2 = prec_words;
    } else {
      k--;
      prec_change = 1;
    }
    mpmul(b, b, sk0, prec_words);
    mpsub(a, sk0, sk1, prec_words);
    prec_words = nw1;
    mpdiv(sk1, b, sk0, prec_words);
    mpmuld(sk0, 0.5, 0, sk1, prec_words);
    prec_words = nw2;
    mpadd(b, sk1, b, prec_words);
    //the above line needs to change if mpadd is not safe for 
    // same variable input/output.

    if(k == mq - nit && iq == 0) {
      iq = 1;
      prec_change = 0;
    }
  } // end for
#if 0
  printf("mpsqrt[4] b=\n");
  cout << b << endl;
#endif
  // Restore original precision level
  prec_words = nws;
  mproun(b);
#if 0
  printf("mpsqrt[5] after mproun b=\n");
  cout << b << endl;
#endif
  if(debug_level >= 7) {
    cerr << "MPSQRT done." << endl;
  } 
}
