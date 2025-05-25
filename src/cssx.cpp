#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpcssx(const mp_real& a, const mp_real& pi,
		     mp_real& x, mp_real& y)
{
  /**
   * This computes the cosine and sine of the MP number A
   * and returns the two MP results in X and Y, respectively.
   *  PI is the MP value of Pi computed by a previous call to MPPI
   * or MPPIX.  Before calling MPCSSX, the arrays mpuu1 and
   * mpuu2 must be initialized by calling MPINIX.  For modest levels
   * of precision, use MPCSSN.  The last word of the result is not
   * reliable.
   *
   * This routine emplys a complex arithmetic version of the
   * scheme found in MPEXPX.
   */
  const double cl2 = 1.4426950408889633, cpi = 3.141592653589793;
  const int nit = 1;
  int prec_words = mp::prec_words;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(x);
    zero(y);
    return;
  }
  if(debug_level >= 5) cerr << "MPCSSX I" << endl;

  int na = std::min(int(std::abs(a[1])), prec_words);
  int ncr = 1 << mpmcrx;

  // Check if precision level is too low to justify advanced routine.

  if(prec_words <= ncr) {
    mpcssn(a, pi, x, y, prec_words); return;
  }

  // Check if input is zero.
  if(!na) {
    //x = 1.0, y = zero.
    x[1] = 1.0; x[2] = 0.0; x[3] = 1.0;
    zero(y);
    return;
  }
  
  // Check if Pi has been precomputed.
  
  double t1;
  int n1;
  mpmdc(pi, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - cpi) > mprx2) {
    if(MPKER[30] != 0) {
      cerr << "*** MPCSSX: PI must be precomputed." << endl;
      error_no = 30;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  int n6 = prec_words+6;
  mp_complex sk0(0.0, 0.0, n6), sk1(0.0, 0.0, n6);
  mp_complex sk3(0.0, 0.0, n6), sk2(0.0, 0.0, n6);
  mp_real f1(1.0, 8);
  int nws = prec_words;
  
  // Reduce argument to betwwen -Pi and Pi.

  mpmuld(pi, 2.0, 0, sk0.real, prec_words);
  mpdivx(a, sk0.real, sk1.real, prec_words);
  mpnint(sk1.real, sk2.real, prec_words);
  mpmulx(sk2.real, sk0.real, sk1.real, prec_words);
  mpsub(a, sk1.real, sk0.real, prec_words);

  // If reduced argument is close to 0, pi/2, -pi/2, pi or -pi, call mpcssn.
  mpmdc (sk0.real, t1, n1, prec_words);
  if (n1 < -2*mpnbt) {
      mpcssn (sk0.real, pi, x, y, prec_words);
      return;
  }
  double t2 = ldexp(t1, n1) / cpi;
  if ( fabs (t2) < 1.0e-10 || fabs (t2 - 0.5) < 1.0e-10 ||
       fabs (t2 + 0.5) < 1.0e-10 || fabs (t2 - 1.0) < 1.0e-10 ||
       fabs (t2 + 1.0) < 1.0e-10 ) {
      mpcssn (sk0.real, pi, x, y, prec_words);
      return;
  }

  // Determine the least integer MQ such that 2 ^ MQ >= prec_words.

  int mq = int(cl2 * log(double(nws)) + 1.0 - mprxx);
  mpeq(f1, sk2.real, prec_words);
  // Compute initial approximation to [Cos (A), Sin (A)].
  
  prec_words = ncr+1;
  mpcssn(sk0.real, pi, sk3.real, sk3.imag, prec_words);
  int iq =0, k, prec_change=1;

  // Perform the Newton-Raphson iteration with a dynamically changing precision
  // level prec_words.
  
  for(k=mpmcrx+1; k<= mq; k++) {
    if(prec_change) {
      prec_words = std::min(2*prec_words-2, nws) + 1;
    } else {
      prec_change = 1;
    }
    mpangx(sk3.real, sk3.imag, pi, sk1.real);
    mpsub(sk0.real, sk1.real, sk2.imag, prec_words);
    mp_complex::mpcmulx(sk3, sk2, sk1, prec_words);
    mp_complex::mpceq(sk1, sk3, prec_words);
    if(k >= mq - nit && iq == 0) {
      iq = 1;
      prec_change = 0;
      k--;
    }
  }
  
  // The final (cos, sin) result must be normalized to magnitude 1.
  mpsqx(sk3.real, sk0.real, prec_words);
  mpsqx(sk3.imag, sk0.imag, prec_words);
  mpadd(sk0.real, sk0.imag, sk1.real, prec_words);
  mpsqrtx(sk1.real, sk2.real, prec_words);
  mpdivx(sk3.real, sk2.real, sk0.real, prec_words);
  mpdivx(sk3.imag, sk2.real, sk0.imag, prec_words);

  prec_words = nws;
  mpeq(sk0.real, x, prec_words);
  mpeq(sk0.imag, y, prec_words);
  
  if (debug_level >= 5) cerr << "MPCSSX 0" << endl;
}

