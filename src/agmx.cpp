#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpagmx(mp_real &a, mp_real& b, int prec_words)
{
  /**
   * This performs the arithmetic-geometric mean (AGM) iterations.
   * this routine is called by MPLOGX.  It is not intended to be 
   * called directly by the user.
   */

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(a); zero(b);
    return;
  }
  int n6 = prec_words + 6;
  mp_real sk0(0.0, n6), sk1(0.0, n6);
  int l1 = 0;
  double s1;
  sk0[2] = 10000000.0;// high value to force second iteration.

  do {
    l1++;
    if (l1 == 50) {
      if(MPKER[5] != 0) {
	cerr <<"*** MPAGMX: Iteration limit exceeded." << endl;
	error_no = 5;
	if(MPKER[error_no] == 2) mpabrt();
      }
      break;
    }
    
    s1 = sk0[2];
    mpadd(a, b, sk0, prec_words);
    mpmuld(sk0, 0.5, 0, sk1, prec_words);
    mpmulx(a, b, sk0, prec_words);
    mpsqrtx(sk0, b, prec_words);
    mpeq(sk1, a, prec_words);
    mpsub(a, b, sk0, prec_words);
    
    // Check for convergence.
  }
  while(sk0[1] != 0.0 && (sk0[2] < s1 || sk0[2] >= -2));
  if (debug_level >= 6) 
    cerr << "MPAGMX: Iteration = " << l1 << ", Tol. Achieved = " << sk0[2] << endl;
  return;
}
 
