#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpsqrtx(const mp_real& a, mp_real& b, int prec_words)
{
  /* This computes the square root of the MP number A and returns the MP
   * result in B.  Before calling MPSQRTX, the arrays mpuu1 and mpuu3 must
   * be initialized by calling MPINIX.  For modest levels of precision,
   * use MPSQRT.
   */
  /**
   * This routine calles MPSQT to obtain an initial
   * approximation.  The newton iteration which
   * converges to 1 / Sqrt(A) used by this routine is :
   *
   *  x_{k+1} = x_k + 0.5 * (1 - A * x_k^2) * x_k
   *
   * Where the multiplication () * x_k is performed at half precision.
   *
   * The finial iteration uses a different iteration
   * and is due to A. Karp : 
   *
   *  x_n = (A * x_{n-1}) + 0.5 * [A - (A * x_{n-1}) ^ 2] * x_{n-1}
   *
   * Where the multiplication (A * x_{n-1}) and [] * x_{n-1} are
   * performed at half precision.
   *
   */
  
  const double cl2 = 1.4426950408889633;
  const int nit = 3;
  
  if (error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 6) {
    cerr << "MPSQRTX I" << endl;
  }
  
  double ia = sign(1.0, a[1]);
  int na = std::min(int(std::abs(a[1])), prec_words);
  int ncr = 1 << mpmcrx;
  
  if(na == 0) {
    zero(b);
    return;
  }
  if(ia < 0.0) {
    if(MPKER[71] != 0) {
      cerr <<"*** MPSQRTX: Argument is negative." << endl;
      error_no = 71;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  // Check if precision level is too low to justify the advanced routine.
  
  if(prec_words < ncr) {
    mpsqrt(a, b, prec_words);
    return;
  }
  
  int nws = prec_words;
  int n6 = prec_words+6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  mp_real f(1.0, 8);
  // Determine the least integer MQ such that 2 ^ MQ > MPNW.

  double t1 = double(prec_words);
  int mq = int(cl2 * log(t1) + 1.0 - mprxx);
  
  // Compute the initial approximation of 1 / Sqrt(A).
  
  prec_words = ncr + 1;
  mpsqrt(a, sk0, prec_words);
  mpdiv(sk0, a, b, prec_words);
  int iq = 0, nw1=0, nw2=0, k;
  int prec_change = 1;
  
  // Perform the Newton-Raphson iteration described above with
  // a dynamically changing precision level MPNW (one greater than powers of two
  for(k=mpmcrx + 1; k < mq; k++) {
    if(prec_change) {
      nw1 = prec_words;
      prec_words = std::min(2*prec_words - 2, nws) +1;
      nw2 = prec_words;
    } else {
      prec_change = 1;
    }
    mpsqx(b, sk0, prec_words);
    mpmulx(a, sk0, sk1, prec_words);
    mpsub(f, sk1, sk0, prec_words);
    prec_words = nw1;
    mpmulx(b, sk0, sk1, prec_words);
    mpmuld(sk1, 0.5, 0, sk0, prec_words);
    prec_words = nw2;
    mpadd(b, sk0, b, prec_words);
    if(k >= mq - nit && iq == 0) {
      iq = 1;
      prec_change = 0;
      k--;
    }
  }
  
  // Perform last iteration using Karp's trick

  mpmulx(a, b, sk0, prec_words);
  nw1 = prec_words;
  prec_words = std::min(2*prec_words-2, nws) + 1;
  nw2 = prec_words;
  mpsqx(sk0, sk1, prec_words);
  mpsub(a, sk1, sk2, prec_words);
  prec_words = nw1;
  mpmulx(sk2, b, sk1, prec_words);
  mpmuld(sk1, 0.5, 0, sk2, prec_words);
  prec_words = nw2;
  mpadd(sk0, sk2, b, prec_words);
 
  //Restore original precision level.
   
  prec_words = nws;
  mproun(b);
}

