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

void mp_real::mplog(const mp_real& a, const mp_real& al2, mp_real& b, 
                    int prec_words, int nit/* = 3*/)
{
  /**
   * This function computest the natural logarithm of the MP number a
   * and puts the result in B.  AL2 is the MP value of log(2), 
   * produced by a prior call to mplog.  For extra high levels of precision, 
   * use mplogx. 
   * 
   * For this routine, the last word of the result is not reliable.
   * debug output starts with debug_level >= 7.
   * 
   * The Taylor series for Log converges more slowly than that of exp.
   * thus this routine does not employ taylor series, but instead computes
   * logarithms by solving exp(b) = a, using the following newton iteration,
   * which converges to log(a) == b.
   * 
   * x_{k+1} = x_k + [a - exp(x_k)] / exp(x_k);
   * 
   * These iterations are performed with a maximum precision level prec_words
   * that is dynamically changed, approximately doubling with each iteration.
   * 
   * See the comment about the int NIT in mpdiv.
   * NIT is an optional argument.  In general it should not be used.
   * The argument is only set (to zero) when creating Log(2) for the 
   * first time. This increases the precision of the computation of Log(2).
   */
  const double alt = 0.693147180559945309;
  const double cl2 = 1.4426950408889633;
  double ia;
  int na, n1, n2;
  double t1, t2;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if (debug_level >= 6) cerr << "MPLOG I" << a << endl;
  
  ia = sign(1.0, a[1]);
  na = std::min(int(std::abs(a[1])), prec_words);
  
  if(ia < 0 || na == 0) { // negative or zero.
    if(MPKER[50] != 0) {
      cerr << "*** MPLOG: Argument is less than or equal to zero." << endl;
      error_no = 50;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  //        Unless the input is close to 2, log(2) must have been
  //        Initialized previously.
  
  mpmdc(a, t1, n1, prec_words);
  if((std::abs(t1 - 2.0) > 0.001)  || (n1 != 0)) {
    mpmdc(al2, t2, n2, prec_words);
    if(n2 != -mpnbt || (std::abs(t2 * mprdx - alt) > mprx2)) {
      if(MPKER[51] != 0) {
        cerr << "*** MPLOG: LOG (2) must be precomputed." << endl;
        error_no = 51;
        if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }
  
  // Check if input is exactly one.
  if((a[1] == 1.0)  && (a[2] == 0.0) && (a[3] == 1.0)) {
    zero(b);
    return;
  }
  
  int n6 = prec_words + 6;
  int nws = prec_words;
  int iq = 0;
  int prec_change = 0, mq, k;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  
  //         Determine the lest integer MQ such that 2 ^ MQ >= prec_words.
  
  t2 = nws;
  mq = int(cl2 * log(t2) + 1.0 -mprxx);
  
  // Compute initial approximation of LOG(A);
  
  t1 = log(t1) + n1*alt;
  mpdmc(t1, 0, b, prec_words);
  prec_words = 3;

  //cerr << "A" << endl;

  //        Perform the Newton-Raphson iteration described above, 
  //        Changing precision level prec_words dynamically. (prec_words =
  //        one greater than powers of two)
  
  for(k=1;k<=mq;k++) {
    if(prec_change)
      prec_words = std::min(2*prec_words - 2, nws)+1;
    else
      prec_change = 1;
    mpexp(b, al2, sk0, prec_words);
    mpsub(a, sk0, sk1, prec_words);
    mpdiv(sk1, sk0, sk2, prec_words);
    mpadd(b, sk2, sk1, prec_words);
    mpeq(sk1, b, prec_words);
    if((k == mq-nit) && !iq) {
      iq = 1;
      k--;
      prec_change = 0;
      if(nit == 0) //used only for computing log(2) from mpinit.
        k-= 3; //This gives additional iterations. 3 was arbitrary.
    }
  } // end for;
  //cerr << "B" << endl;
  
  //        Restore original precision level.
  
  prec_words = nws;
  //cerr << "C" << endl;
  //cerr << "sk1[0] = " << sk1[0] << endl;
  //cerr << "b[0] = " << b[0] << endl;
  //cerr << "sk1[1] = " << sk1[1] << endl;
  //cerr << "b[1] = " << b[1] << endl;
  if (sk1[1] > b[1])
    mproun(b);
  //cerr << "D" << endl;
  return;
} // End MPLOG.

