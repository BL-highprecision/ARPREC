#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpcshx(const mp_real& a, const mp_real& pi, 
		     const mp_real& al2, mp_real& x, mp_real& y)
{
  /** 
   * This computes the hyperbolic cosine and sine of the MP
   * number A and returns the two MP results in X and Y,
   * respectively.  PI is the MP value of Pi computed by a previous
   * call to MPPI or MPPIX.   AL2 is the MP value of Log(2),
   * computed by a previous call to MPLOG or MPLOGX.  Before calling
   *  MPCSHX, the arrays mpuu1 and mpuu2 must be initialized by 
   * calling MPINIX.  For modest levels of precision, use MPCSSH.
   */
  int prec_words = mp::prec_words;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(x);
    zero(y);
    return;
  }
  if(debug_level >= 5) cerr << "MPCSHX I" << endl;

  int n6 = prec_words + 6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6), f(1.0, 8);

  mpexpx(a, pi, al2, sk0);
  mpdivx(f, sk0, sk1, prec_words);
  mpadd(sk0, sk1, sk2, prec_words);
  mpmuld(sk2, 0.5, 0, x, prec_words);
  mpsub(sk0, sk1, sk2, prec_words);
  mpmuld(sk2, 0.5, 0, y, prec_words);
  
  if(debug_level >= 5) cerr << "MPCSHX 0" << endl;
  return;
}

