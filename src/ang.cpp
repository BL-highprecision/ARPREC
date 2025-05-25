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

void mp_real::mpang(const mp_real& x, const mp_real& y, 
                    const mp_real& pi, mp_real& a)
{
  /**
   * This computes the MP angle A subtended by the MP pair (x, y)
   * considered as a point in the x-y plane.  This is more useful than
   * an arctan or arcsin routine, since it places the result correctly
   * in the full circle, i.e. -Pi < A <= Pi.  PI is the MP value of
   * Pi computed by a precious call to MPPI.  For extra high levels of
   * precision, use MPANGX.  The last word of the result is not reliable.
   * Debug output starts with debug_level == 5.
   *
   * Space required for A : prec_words + 4 (double) cells.
   *
   * The Taylor series for ArcSin converges much more slowly than that
   * of Sin.  Thus this routine does not empoly Taylor series, but instead
   * computes Arccos or Arcsin by solving Cos(a) = x or Sin(a) = y 
   * using one of the following Newton iterations, both of which
   * converge to a :
   * 
   *  z_{k=1} = z_k - [x - Cos (z_k)] / Sin (z_k)
   *  z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)
   *
   * The first is selected if Abs (x) <= Abs (y); otherwise the second is 
   * used.  These iterations are performed with a maximum precision level
   * prec_words that is dynamically changed, approximately doubling with each 
   * iteration.
   * See the comment about the parameter NIT in MPDIVX.
   */
  const double cl2 = 1.4426950408889633, cpi = 3.141592653589793;
  const int nit = 3;
  double ix, iy, t1, t2;
  int nx, ny, n1, n2;
  int prec_words = mp::prec_words;
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(a);
    return;
  }
  if (debug_level >= 5)
    cerr << "MPANG: x = " << x << ", y = " << y << endl;

  ix = sign(1.0, x[1]);
  nx = std::min(int(std::abs(x[1])), prec_words);
  iy = sign(1.0, y[1]);
  ny = std::min(int(std::abs(y[1])), prec_words);
  
  //  Check if both x and y are zero.
  if(!nx && !ny) {
    if(MPKER[7] != 0) {
      cerr <<"*** MPANG : Both arguments are zero." << endl;
      error_no = 7;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  // Check if Pi has been precomputed
  
  mpmdc(pi, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - cpi) > mprx2) {
    if(MPKER[8] != 0) {
      cerr << "*** MPANG: Pi must be precomputed." << endl;
      error_no = 8;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  // Check if one of x or y is zero
  if(!nx) {
    if (iy > 0.0) 
      mpmuld(pi, 0.5, 0, a, prec_words);
    else
      mpmuld(pi, -0.5, 0, a, prec_words);
    return;
  } else if(!ny) {
    if (ix > 0.0) 
      zero(a);
    else
      mpeq(pi, a, prec_words);
    return;
  }
  
  int nws = prec_words;
  prec_words++;
  int n5 = prec_words+5;
  mp_real sk0(0.0, n5), sk1(0.0, n5), sk2(0.0, n5), sk3(0.0, n5), sk4(0.0, n5);
  int mq, k, kk;
  int iq, prec_change;
  double t3;

  // Determine the least integer mq such that 2^mq >= prec_words
  
  t1 = double(nws);
  mq = int(cl2*log(t1) + 1.0 -mprxx);
  
  // Normalize x and y so that x^2 + y^2 == 1.

  mpmul(x, x, sk0, prec_words);
  mpmul(y, y, sk1, prec_words);
  mpadd(sk0, sk1, sk2, prec_words);
  mpsqrt(sk2, sk3, prec_words);
  mpdiv(x, sk3, sk1, prec_words); // sk1 holds scaled x
  mpdiv(y, sk3, sk2, prec_words); // sk2 holds scaled y

  // Compute initial approximation of the angle.
  mpmdc(sk1, t1, n1, prec_words);
  mpmdc(sk2, t2, n2, prec_words);
  n1 = std::max(n1, -66);
  n2 = std::max(n2, -66);
  t1 *= pow(2.0, double(n1));
  t2 *= pow(2.0, double(n2));
  t3 = atan2(t2, t1);
  mpdmc(t3, 0, a, prec_words);
  
  // The smaller of x or y will be used from now on to measure convergence.
  //This selects the Newton iteration of the two listed above that has
  // The largest denominator.
  
  if (std::abs(t1) <= std::abs(t2)) {
    kk = 1;
    mpeq(sk1, sk0, prec_words);
  } else {
    kk = 2;
    mpeq(sk2, sk0, prec_words);
  }
  
  prec_words = 3;
  iq = 0;
  prec_change = 0; 
  
  // Perform the Newton-Raphson iteration described above with a dynamically
  // changing precision level prec_words (one greater than powers of 2).

  for(k=1;k<=mq;k++) {
    if(prec_change) 
      prec_words = std::min(2*prec_words-2, nws) + 1;
    else
      prec_change = 1;
    //most work done here.
    mpcssn(a, pi, sk1, sk2, prec_words);

    if(kk == 1) {
      mpsub(sk0, sk1, sk3, prec_words);
      mpdiv(sk3, sk2, sk4, prec_words);
      mpsub(a, sk4, sk1, prec_words);
    } else {
      mpsub(sk0, sk2, sk3, prec_words);
      mpdiv(sk3, sk1, sk4, prec_words);
      mpadd(a, sk4, sk1, prec_words);
    }
    mpeq(sk1, a, prec_words);
    if(k == mq - nit && !iq) {
      iq = 1; 
      prec_change = 0;
      k--;
    }
  }
  
  // restore original precision level.
  prec_words = nws;
  mproun(a);
  
  if (debug_level >= 5) cerr << "MPANG: done :" << a << endl;
  return;
}

