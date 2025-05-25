/*
 * src/mpreal4.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Additional routines, mostly those for which speed is less important.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mppi(mp_real& pi)
{
  /**
   * This computes Pi to available precision (MPNW mantissa words). For
   * extra high levels of precision, use MPPIX.  The last word of the
   * result is not reliable.  Debug output starts with debug_level = 7.
   *
   * Required space in pi : prec_words + 4 cells.
   * 
   * The algorithm that is used for computing Pi, which is due to Salamin
   * and Brent, is as follows: 
   *
   * Set A_0 = 1, B_0 = 1/Sqrt(2), and D_0 = Sqrt(2) - 1/2.
   * 
   * Then from k = 1 iteratoe the following operations:
   * 
   * A_k = 0.5 * (A_{k-1} + B_{k-1})
   * B_l = Sqrt (A_{k-1} * b_{k-1})
   * D_k = D_{k-1} - 2^k * (A_k - B_k) ^2.
   * 
   * The P_k = (A_k + B_k) ^2 / D_k  converges quadratically to Pi.
   * In other words, each iteration approximately doubles the 
   * number of correct digits, providing all iterations are
   * done with the maximum precision.
   */
  const double cl2 = 1.4426950408889633;
  int prec_words = mp::prec_words;
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(pi);
    return;
  }

  // Perform calculations to one extra word accuracy.
  double t1;
  int mq, k;
  int nws = prec_words;
  int n6 = prec_words +6;
  prec_words++;
  mp_real sk0(1.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  mp_real sk3(0.0, n6), sk4(0.0, n6), f(2.0, 8);
    
  // Determine the number of iterations required for the
  // given precision level.
  
  t1 = nws * log10(mpbdx);
  mq = int(cl2 * (log(t1) - 1.0) + 1.0);
#if 0
  printf("mppi: prec_words %d, mpbdx %22.18e, mq %d\n", prec_words, mpbdx, mq);
#endif 
  // Initialize as above.
  mpsqrt(f, sk2, prec_words);
#if 0
  printf("mppi[1] after mpsqrt: sk2, f\n"); 
  cout << sk2 << endl;
  cout << f << endl;
#endif
  mpmuld(sk2, 0.5, 0, sk1, prec_words);
  // f = 1/2;
  f[2] = -1; f[3] = 0.5 * mpbdx;
  mpsub(sk2, f, sk4, prec_words);
#if 0
  printf("mppi[2]: sk2, f, sk4 \n"); 
  cout << sk2 << endl;
  cout << f << endl;
  cout << sk4 << endl;
#endif  
  // Perform iterations as described above.
  
  for(k=1; k<=mq ; k++) {
    mpadd(sk0, sk1, sk2, prec_words);
    mpmul(sk0, sk1, sk3, prec_words);
    mpsqrt(sk3, sk1, prec_words);
    mpmuld(sk2, 0.5, 0, sk0, prec_words);
    mpsub(sk0, sk1, sk2, prec_words);
    mpmul(sk2, sk2, sk3, prec_words);
    t1 = ldexp(1.0, k);
    mpmuld(sk3, t1, 0, sk2, prec_words);
    mpsub(sk4, sk2, sk3, prec_words);
    mpeq(sk3, sk4, prec_words);
  }
#if 0
  printf("mppi[3]: sk4 \n");
  cout << sk4 << endl;
#endif
  // Complete computation.
  
  mpadd(sk0, sk1, sk2, prec_words);
  mpmul(sk2, sk2, sk3, prec_words);
  mpdiv(sk3, sk4, sk2, prec_words);
  mpeq(sk2, pi, prec_words);
  
  // Restore original precision level.
  prec_words = nws;
  mproun(pi);
  
  if(debug_level >= 7) cerr <<"Computed Pi: "<< pi << endl;
}

