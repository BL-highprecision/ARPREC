#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpangx(const mp_real& x, const mp_real& y, 
		     const mp_real& pi, mp_real& a)
{
  /**
   * This computes the MP angle A subtended by the MP pair (X, Y)
   * considered as a point in the x-y plane.  This is more usefull than
   * an arctan or arcsin routine, since it places the result correctly in
   * the full circle, i.e. -Pi < A <= Pi.  PI is the MP value of Pi computed
   * by a previous call to MPPI or MPPIX.  Before calling MPANGX, the arrays
   * mpuu1 and mpuu2 must be initialized by calling MPINIX.  For modest
   * levels of precision, use MPANG.  The last word of the result
   * is not reliable.
   *
   * This routine employs a complex arithmetic version of the MPLOGX algorithm.
   */

  const double cpi = 3.141592653589793;
  int prec_words = mp::prec_words;

  if(error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(a);
    return;
  }
  if (debug_level >= 6) cerr << "MPANGX I" << endl;
  
  int ix = sign(1, int(x[1]));
  int nx = std::min(int(std::abs(x[1])), prec_words);
  int iy = sign(1, int(y[1]));
  int ny = std::min(int(std::abs(y[1])), prec_words);
  int ncr = 1 << mpmcrx;

  // Check if precision level it too low to justify the advanced routine.
  
  if(prec_words <= ncr) {
    mpang(x, y, pi, a); return;
  }

  // Check if both X and Y are zero.
  if(!nx && !ny) {
    if(MPKER[9] != 0) {
      cerr << "*** MPANGX: Both arguments are zero." << endl;
      error_no = 9;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  // Check if Pi has been precomputed.
  
  double t1;
  int n1;
  mpmdc(pi, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - cpi) > mprx2) {
    if(MPKER[10] != 0) {
      cerr << "*** MPANGX: PI must be precomputed." << endl;
      error_no = 10;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  // Check if one of X or Y is zero.
  if(nx == 0) {
    if(iy > 0) 
      mpmuld(pi, 0.5, 0, a, prec_words);
    else
      mpmuld(pi, -0.5, 0, a, prec_words);
    return;
  } else if(ny == 0) {
    if(ix > 0) 
      zero(a);
    else
      mpeq(pi, a, prec_words);
    return;
  }
  
  int n6 = prec_words+6;
  mp_complex sk0(0.0, 0.0, n6), sk1(0.0, 0.0, n6), sk2(0.0, 0.0, n6);
  mp_complex f1(1.0, 0.0, 8), f4(4.0, 0.0, 8);
  zero(f1.imag);
  zero(f4.imag);

  // Multiply the input by a large power of two.

  mpmdc(x, t1, n1, prec_words);
  int n2 = mpnbt * (prec_words / 2 + 2) - n1;
  mpmuld(x, 1.0, n2, sk0.real, prec_words);
  mpmuld(y, 1.0, n2, sk0.imag, prec_words);

  // Perform AGM iterations.

  mp_complex::mpceq(f1, sk1, prec_words);
  mp_complex::mpcdivx(f4, sk0, sk2, prec_words);
  mp_complex::mpcagx(sk1, sk2);

  // Compute A = Imag (Pi / (2*Z)), where Z is the limit of the complex
  // AGM.

  mp_complex::mpcmuld(sk1, 2.0, 0, sk0, prec_words);
  mpeq(pi, sk2.real, prec_words);
  zero(sk2.imag);
  mp_complex::mpcdivx(sk2, sk0, sk1, prec_words);
  mpeq(sk1.imag, a, prec_words);
  
  return;
}

