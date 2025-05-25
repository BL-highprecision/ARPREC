#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpsqx(const mp_real& a, mp_real& b, int prec_words)
{
  /**
   * This routine squares the MP number A to yield the MP product B. 
   * Before calling MPSQX, the arrays mpuu1 and mpuu3 must be initialized
   * by calling mpinix.  For modest levels of precision, use MPMUL.
   *
   * This function uses the same FFT technique as MPMULX.  It is 
   * faster because only one forward FFT has to be computed.
   */
  

  double t1, t2;
  int ia, na, ncr, i, i2=1, nn, nb, nx;

  if(error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if (debug_level >= 8) cerr << "MPSQX I" << endl;
  
  ia = sign(1, int(a[1]));
  na = std::min(int(std::abs(a[1])), prec_words);
  ncr = 1 << mpmcrx; // pow(2.0, mpmcrx);
  if(!na) {
    // one of the inputs is zero -- result is zero.
    zero(b);
    return;
  }

  // Check if precision level of one of the arguments is 
  // enough to justify the advanced routine.
  
  if(na <= ncr) {
    mpmul(a, a, b, prec_words);
    return;
  }

  //now allocate arrays.
  double *d1 = new double[4*prec_words+8];
  double *d3 = new double[8*prec_words+16];
  //mp_down14 represents value to multiply by
  //to get one fourth of 48 bits below the decimal point.
  const double mp_down14 = 1.0 / 4096.0;
  const double mp_down24 = mprbx;
  const double mp_down34 = mp_down14 * mp_down24;
  const double mp_up14 = 4096.0;
  const double mp_up24 = mpbbx;
  const double mp_up34 = mp_up24 * mp_up14;


  // Place the input data in A into the scratch array DD1. 
  // this code alsp splits the input data into fourth sized words.
  // (mpnbt/4) bits for each new word (maximum)
  
  for(i=0;i<na;i++) {
    i2 = 4*i;
    t1 = a[i+FST_M];
    t2 = FLOOR_POSITIVE(mp_down34*t1);
    d1[i2] = t2;
    t1 -= mp_up34 * t2;
    t2 = int(mp_down24*t1);
    d1[i2+1] = t2;
    t1 -= mp_up24 * t2;
    t2 = int(mp_down14*t1);
    d1[i2+2] = t2;
    t1 -= mp_up14 * t2;
    d1[i2+3] = t1;
  }

  // Call the convolution

  nn = 4 * na;
  nx = int(ANINT(sqrt(3.0 * nn) + mprxx));
  mplconv(1, nn, nx, d1, d1, d3);

  // Recombine words and release carries.

  nb = std::min(na + na, prec_words+3);
  int nb1 = nb - 1;
  d1[1] = nb; //result always positive
  d1[2] = a[2] + a[2] + 1;

  d1[3] = d3[0] * mp_up24 + d3[1] * mp_up14 + d3[2];
  d1[nb+FST_M+1] = d1[nb+FST_M+2] = 0.0;
  
  for(i=0;i<nb1;i++) {
    i2 = i * 4 + 3;
    t1 = d3[i2];
    t2 = FLOOR_POSITIVE(t1 * mp_down14);
    t1 -= t2 * mp_up14;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] = t1 * mp_up34;
    t1 = d3[i2+1];
    t2 = int(t1 * mp_down24);
    t1 -= t2 * mp_up24;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] += t1 * mp_up24;
    t1 = d3[i2+2];
    t2 = int(t1 * mp_down34);
    t1 -= t2 * mp_up34;
    d1[i+FST_M] += t2;
    d1[i+FST_M+1] += t1 * mp_up14 + d3[i2+3];
  }
  //one more to do.
  t1 = d3[i2+4];
  t2 = FLOOR_POSITIVE(t1 * mp_down14);
  t1 -= t2 * mp_up14;
  d1[i+FST_M] += t2;
  d1[i+FST_M+1] = t1 * mp_up34;
  
  //now go on.
  int d_add=0;
  i=0;
  //eliminate leading zeros (except where something
  // will probably carry in)
  while(i<(nb1-3) && d1[i+FST_M] == 0.0 && d1[i+FST_M+1] < (mpbdx-2.0)) {
    i++;
  }
  if(i) {
    d1[2] -= double(i);
    d1[1] = sign(1.0, d1[1]) * (std::abs(d1[1])  - double(i));
    d1[i+2] = d1[2];
    d1[i+1] = d1[1];
    d1 += i;
    d_add -= i;
  }

  // Fix up the result;
  mpnorm(d1, b, prec_words);
  
  if(debug_level >= 8) cerr << "MPSQX 0" << endl;
  delete [] (d1+d_add);
  delete [] d3;
}

