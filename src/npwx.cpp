#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpnpwx(const mp_real& a, int n, mp_real& b, int prec_words)
{
  /**
   * This computes the N-th power of the MP number A and places
   * the MP result in B.  When N is zero, 1 is returned.  When N is negative,
   * the reciprocal of A ^ |N| is returned.  Before calling MPNPWX, the
   * arrays of mpuu1 and uu2 must be initialized by calling mpinix.
   * For modest levels of precision, use MPNWPR.
   *
   *
   * This routine emplys the normal binary method for exponentiation.
   */
  
  const double cl2 = 1.4426950408889633;// 1/log(2)

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 6) {
    cerr << "MPNPWX I" << endl;
  }

  int ncr = 1 << mpmcrx;
  int na = std::min(int(std::abs(a[1])), prec_words);

  // Check if precision level of A is too low to justify the advanced routine.
  
  if(na < ncr && n >= 0 && n <= 4) {
    mpnpwr(a, n, b, prec_words);
    return;
  }
  if(na == 0) {
    if(n >= 0) {
      zero(b);
      return;
    } else {
      // error : 0^N, N <= 0.
      if(MPKER[58] != 0) {
	cerr <<"*** MPNPWX: argument is zero and N is negative or zero." << endl;
	error_no = 58;
	if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }

  int n5 = prec_words+5;
  int nn = std::abs(n);
  mp_real sk0(0.0, n5), sk1(0.0, n5), f1(1.0, 8);

  if(nn == 0) {
    mpeq(f1, b, prec_words);
    return;
  } else if(nn == 1) {
    mpeq(a, b, prec_words);
  } else if(nn == 2) {
    mpsqx(a, b, prec_words);
  } else {
    // full binary exponentiation.
    // Determine the least integer mn such that 2 & mn > nn.
    int mn = int(ANINT(cl2 * log(double(nn)) + 1.0 + mprxx));
    mpeq(a, sk0, prec_words);
    mpeq(f1, b, prec_words);
    int kn = nn, j, kk;

    for(j=1; j <= mn; j++) {
      kk = kn / 2;
      if(kn != 2 * kk) {
	mpmulx(b, sk0, b, prec_words);
      }
      kn = kk;
      if(kn && j < mn)
	mpsqx(sk0, sk0, prec_words);
    }
  }
  
  // Compute the reciprocal if N is negative.
  if(n < 0) {
    mpdivx(f1, b, sk0, prec_words);
    mpeq(sk0, b, prec_words);
  }
  return;
}

