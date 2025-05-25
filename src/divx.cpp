#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpdivx(const mp_real& a, const mp_real& b, mp_real& c, int prec_words)
{
  /**
   * This divides the MP number A by the MP number B and returns the MP
   * result in C.  Before calling MPDIVX, the arrays mpuu1, and mpuu2 must
   * be initizlied by calling mpinix.  for modest levels of precision, use
   * mpdiv.  debug output starts with debug_level = 7.

   * This subroutine emplys the following Newton-Raphson iteration, which
   * converges to 1 / B:
   *
   *  X_{k+1} = X_k + (1 - X_k * B) * X_k
   *
   * where the multiplication () * X_k is performed with only half of 
   * the normal level of precision.  These iterations are performed with
   * a maximum precision level MPNW that is dynamically changing with
   * each iteration.  The final iteration is performed as follows (this
   * is due to A. Karp):
   
   * A / B = (A * X_n) + [A - (A * X_n) * B] * X_n (approx.).

   * where the multiplications A * X_n and [] * X_n are performed with
   * only half the final level of precision.
   *
   * One difficulty iwht this procedure is that errors often accumulate
   * in the trailing mantissa words. This error can be controlled by repeating
   * one of the iterations.  The iteration that is repeated is controlled by
   * setting the parameter NIT below : if NIT = 0, the last iteration is
   * repeated (this is the most effective but the most expensive).  If NIT 
   * is 1, then the next-to-last iteration is repeated, etc.
   *
   *
   *  
   */

  const double cl2 = 1.442695040889633;
  const int nit = 3;
  const int ncr = 1 << (mpmcrx+2);

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if(debug_level >= 7) {
    cerr << "MPDIVX I" << endl;
  }
  int na = std::min(int(std::abs(a[1])), prec_words);
  int nb = std::min(int(std::abs(b[1])), prec_words);
  

  // Check if precision level of divisor is too low
  // to justify the advanced routine.
  
  if(nb <= ncr) {
    mpdiv(a, b, c, prec_words);
    return;
  }


  //check if divisor is zero.
  if(nb == 0) {
    if(MPKER[33] != 0) {
      cerr << "*** MPDIVX: Divisor is zero." << endl;
      error_no = 33;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  // check if dividend is zero.
  if(na == 0) {
    zero(c);
    return;
  }


  int n6 = prec_words + 6;
  int nws = prec_words;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6), f(1.0, 8);
  mp_real c2(0.0, static_cast<int>(c[0])); //allocates c[0] size array for c2

  // Determine the least integer MQ such that 2^MQ >= prec_words.

  int mq = int(cl2 * log(double(prec_words)) + 1.0 - mprxx);

  // Compute the initial approximation of 1 /B to a precision of ncr words.
  prec_words = ncr + 1;
  mpdiv(f, b, c2, prec_words);
  int iq = 0, prec_change = 1, k;
  int nw1=0, nw2=0;
  
  // Perform the Newton-Raphson iterations described above.

  for(k=mpmcrx+2;k<mq; k++) {
    if(prec_change) {
      nw1 = prec_words;
      prec_words = std::min(2*prec_words - 2, nws) + 1;
      nw2 = prec_words;
    } else {
      prec_change = 1;
    }
    mpmulx(b, c2, sk0, prec_words);
    mpsub(f, sk0, sk1, prec_words);
    prec_words = nw1;
    mpmulx(c2, sk1, sk0, prec_words);
    prec_words = nw2;
    mpadd(c2, sk0, c2, prec_words);
    if(k >= mq - nit && !iq) {
      iq = 1;
      prec_change = 0;
      k--;
    }
  }
  
  //Perform last iteration using Karp's trick.
  mpmulx(a, c2, sk0, prec_words);
  nw1 = prec_words;
  prec_words = std::min(2*prec_words - 2, nws) + 1;
  nw2 = prec_words;
  mpmulx(sk0, b, sk1, prec_words);
  mpsub(a, sk1, sk2, prec_words);
  prec_words = nw1;
  mpmulx(sk2, c2, sk1, prec_words);
  prec_words = nw2;
  mpadd(sk0, sk1, c, prec_words); //note: not c2, but c (avoid copy)
  
  // Restore original precision level.
  prec_words = nws;
  mproun(c);

  if(debug_level >= 7) {
    cerr << "MPDIVIX 0" << endl;
  }
  return;
}

