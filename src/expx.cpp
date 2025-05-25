#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpexpx(const mp_real& a, const mp_real& pi, 
		     const mp_real& al2, mp_real& b)
{
  /**
   * This computes the exponential function of the MP number A and 
   * places the MP result in B.  Pi is the MP value of Mp produced
   * by mppix.  AL2 is the MP value of Log(2) produced by a prior call to 
   * MPLOG or MPLOGX.  Before calling MPEXPX, the arrays mpuu1
   * and mpuu2 must be initialized by calling mpinix.  For modest
   * levels of pecision, use MPEXP.
   * The last word of the result is not reliable.  
   *
   * This routine uses the Newton iteration 
   *  
   *  b_{k+1} = b_k [a - log b_k] + b_k
   *
   * with a dynamically changing level of precision. Logs are performed
   * using MPLOGX.  See the comment about the parameter NIT in MPDIVX.
   * The multiplication b_k * [] is performed at roughly half precision.
   *
   */
  const double alt = 0.693147180559945309;
  const double cl2 = 1.4426950408889633,
    cpi = 3.141592653589793238;
  int nit = 3;
  double t1, t2;
  int n1, n2;
  int prec_words = mp::prec_words;

  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 5) cerr << "MPEXPX I" << endl;
  
  int ncr = 1 << (mpmcrx+1);

  // Check if precision level is low to justify the advanced routine.
  
  if (prec_words <= ncr) {
    mpexp(a, al2, b, prec_words);
    return;
  }

  mpmdc(a, t1, n1, prec_words);
  if(n1)
    t1 = ldexp(t1, n1);
  
  // Check if Log(2) has been precomputed.
  
  mpmdc(al2, t2, n2, prec_words);
  if(n2 != -mpnbt || std::abs(t2 * mp::mprdx - alt) > mprx2) {
    if (MPKER[37] != 0) {
      cerr << "*** MPEXPX: LOG(2) must be precomputed." << endl;
      error_no = 37;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  // check if Pi has been precomputed.
  mpmdc(pi, t2, n2, prec_words);
  if(n2 != 0 || fabs(t2 - cpi) > mprx2) {
    if(MPKER[38] != 0) {
      cerr <<"*** MPEXPX: Pi must be precomputed." << endl;
      error_no = 38;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  // Check for overflows and underflows.
  // The constant is roughly (2^26 * 48)* log(2),
  // Which allows a result with (word) exponent up to 2^26.
  if(fabs(t1) > 2232783353.0) {
    if(t1 > 0.0) {
      if(MPKER[39] != 0) {
	cerr << "*** MPEXPX: Argument is too large" << endl;
	error_no = 39;
	if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    } else {
      zero(b);
      return;
    }
  }

  int n6 = prec_words + 6;
  int nws = prec_words;
  mp_real sk1(0.0, n6), sk2(0.0, n6);

  // Determine the least integer MQ such that 2 ^ MQ >= prec_words.
  
  int mq = int(cl2 * log(double(nws)) + 1.0 - mprxx);

  // compute initial approximation to exp(a).
  
  prec_words = ncr+1;
  mpexp(a, al2, b, prec_words);
  int k, prec_change = 1, iq = 0, nws1, nws2;
  //step up precision
  for(k=mpmcrx+2;k<=mq;k++) {
    nws1 = std::min(prec_words+5, nws) +1;
    if(prec_change) {
      prec_words = std::min(2 * prec_words -2, nws) + 1;
    } else {
      prec_change = 1;
    }
    nws2 = prec_words;
    mplogx(b, pi, al2, sk1, prec_words);
    mpsub(a, sk1, sk2, prec_words);
    prec_words = nws1;
    mpmulx(b, sk2, sk2, prec_words);
    prec_words = nws2;
    //step up precision before add.
    prec_words = std::min(prec_words + 5, nws) + 1;
    mpadd(b, sk2, b, prec_words);
    prec_words = nws2;
    if(!iq && mq -nit >= k) {
      k--;
      prec_change = 0;
      iq = 1;
    }
  }
  if(debug_level >= 6) cerr << "MPEXPX 0" << endl;
  prec_words = nws;
  mproun(b);
}

