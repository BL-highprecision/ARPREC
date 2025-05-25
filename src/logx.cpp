#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mplogx(const mp_real& a, const mp_real& pi, const mp_real& al2, 
		     mp_real& b, int prec_words)
{
  /**
   * This computes the natrural logarithm of the MP number a and places
   * the MP result in B.  PI is the MP value of Pi produced previously
   * by a call to MPPI or MPPIX.  AL2 is the MP value of Log(2), produced by a
   * previous call to MPLOG or MPLOGX.  Before calling MPLOGX, the
   * mpuu1 and mpuu2, and mpuu3 arrays must be initialized by calling
   * MPINIX. For modest levels of precision, use mplog.  the last 
   * word of the result is not reliable.  Debug tarts with debug_level = 6.
   *
   * This routine uses the following algorithm, which is due to Salamin.
   * If a is extremeley close to 1, use a Taylor series.  Otherwise,
   * Select n such that z = x 2^n is at least 2^m, where m is the number
   * of bits of desired precision in the result.  Then,
   *
   *   Log(x) = Pi  / [2 AGM (1, 4/x) ].
   */
  double t1, t2, tn;
  const double alt = 0.693147180559945309, cpi = 3.141592653589793;
  const int mzl = -5;
  int it2, n1;
  
  if(error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if (debug_level >= 6) cerr << "MPLOGX I" << endl;

  int ia = sign(1, int(a[1]));
  int na = std::min(std::abs(int(a[1])), prec_words);
  int ncr = 1 << (mpmcrx-2); //This version of log is faster at
			// smaller number of words.
  int n2;

  // Check if precision level is too low to justify the advanced routine.

  if (prec_words <= ncr) {
    mplog(a, al2, b, prec_words);
    return;
  }
  
  if(ia < 0 || na == 0) {
    //input is less than or equal to zero.
    if(MPKER[52] != 0) {
      cerr <<"*** MPLOGX: Argument is less than or equal to zero." << endl;
      error_no = 52;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  // check if Pi has been precomputed.
  mpmdc(pi, t1, n1, prec_words);
#if 0
  cout << "mplogx[1] After mpmdc: " << " n1 " <<  n1 << endl;
  printf("t1 %22.18e\n", t1);
  cout << "prec_words " << prec_words << endl;
  printf("mprx2 %e, std::abs(t1-cpi) %22.18e\n", mprx2, std::abs(t1-cpi));
#endif
  if(n1 != 0 || std::abs(t1 - cpi) > mprx2) {
    if(MPKER[53] != 0) {
      cerr << "*** MPLOGX: PI must be precomputed." << endl;
      error_no = 53;
      if(MPKER[error_no] == 2) mpabrt();      
    }
    return;
  }

  // Unless the input is 2, Log(2) must have been precomputed.
 
  if(a[1] != 1.0 || a[2] != 0.0 || a[3] != 2.0) {
    it2 = 0;
    mpmdc(al2, t2, n2, prec_words);
    if(n2 != -mpnbt || std::abs(t2 * mprdx - alt) > mprx2) {
      if(MPKER[54] != 0) {
	cerr << "*** MPLOGX: Log (2) must be precomputed." << endl;
	error_no = 54;
      if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  } else {
    it2 = 1;
  }
  
  int nws = prec_words;
  prec_words++;

  int n6 = prec_words + 6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  mp_real sk3(0.0, n6), f1(1.0, 8), f4(4.0, 8);

  // If argument is 1 the result is zero.  If the argement is
  //extremeley close to 1, employ  a Taylor's series instead.

  mpsub(a, f1, sk0, prec_words);
  if(sk0[1] == 0.0) {
    zero(b);
    prec_words = nws;
    return;
  } else if(sk0[2] < mzl) {
    mp_real sk4(0.0, n6);
    mpeq(sk0, sk1, prec_words);
    mpeq(sk1, sk2, prec_words);
    int i1 = 1;
    int tl = int(sk0[2] - prec_words - 1);
    double st, is = 1.0;
    if(debug_level >= 6) cerr <<"Using Taylor series in MPLOGX." << endl;
    do {
      i1++;
      is = -is;
      st  = is * i1;
      mpmulx(sk1, sk2, sk3, prec_words);
      mpeq(sk3, sk2, prec_words);
      mpdivd(sk3, st, 0, sk4, prec_words);
      mpadd(sk0, sk4, sk3, prec_words);
      mpeq(sk3, sk0, prec_words);
    } while(sk2[2] >= tl);

    prec_words = nws;
    mpeq(sk0, b, prec_words);
    return;
  }
  
  // If input is exactly 2, set the exponent to a large value. Otherwise,
  // multiply the input by a large power of two.

  mpmdc(a, t1, n1, prec_words);
  n2 = mpnbt * (prec_words / 2 + 2) - n1;
  tn = n2;
  if(it2 == 1) 
    mpdmc(1.0, n2, sk0, prec_words);
  else 
    mpmuld(a, 1.0, n2, sk0, prec_words);

  // Perform AGM iterations.
  mpeq(f1, sk1, prec_words);
  mpdivx(f4, sk0, sk2, prec_words);
  mpagmx(sk1, sk2, prec_words);
  
  // Compute B = Pi / (2 * A), where A is the limit of the AGM iterations.

  mpmuld(sk1, 2.0, 0, sk0, prec_words);
  mpdivx(pi, sk0, sk1, prec_words);
  // If the input was exactly 2, divide by TN.  Otherwise,
  // subtract TN * Log(2).
  
  if(it2 == 1) {
    mpdivd(sk1, tn, 0, b, prec_words);
  } else {
    mpmuld(al2, tn, 0, sk2, prec_words);
    mpsub(sk1, sk2, b, prec_words);
  }
  prec_words = nws;

  if(debug_level >= 6) cerr << "MPLOGX 0" << endl;
}

