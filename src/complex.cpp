/*
 * src/mpcomplex.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include <arprec/mp_complex.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_complex::mpcmuld(const mp_complex& a, double db, int n, 
			 mp_complex& b, int prec_words)
{
  mp_real::mpmuld(a.real, db, n, b.real, prec_words);
  mp_real::mpmuld(a.imag, db, n, b.imag, prec_words);  
}

void mp_complex::mpcadd(const mp_complex& a, const mp_complex& b,
			mp_complex& c, int prec_words)
{
  mp_real::mpadd(a.real, b.real, c.real, prec_words);
  mp_real::mpadd(a.imag, b.imag, c.imag, prec_words);
}

void mp_complex::mpcsub(const mp_complex& a, const mp_complex& b,
			mp_complex& c, int prec_words)
{
  mp_real::mpsub(a.real, b.real, c.real, prec_words);
  mp_real::mpsub(a.imag, b.imag, c.imag, prec_words);
  return;
}

void mp_complex::mpcmul(const mp_complex& a, const mp_complex& b,
			mp_complex& c, int prec_words)
{
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if(debug_level >= 7)
    cerr << "MPCMUL I" << endl;

  // use 3 multiply method for complex multiplication.
  // This method requires more additions, but with this library,
  // multiplications are more expensive.
  mp_real temp1, temp2, temp4;
  mp_real::mpmul(a.real, b.real, temp1, prec_words);
  mp_real::mpmul(a.imag, b.imag, temp2, prec_words);
  mp_real::mpsub(temp1, temp2, temp4, prec_words);
  mp_real::mpadd(temp2, temp1, temp2, prec_words);
  mp_real temp3;
  mp_real::mpadd(a.real, a.imag, temp1, prec_words);
  mp_real::mpadd(b.real, b.imag, temp3, prec_words);
  mp_real::mpmul(temp1, temp3, temp3, prec_words);
  mp_real::mpsub(temp3, temp2, c.imag, prec_words);
  mp_real::mpeq(temp4, c.real, prec_words);
  return;
}

void mp_complex::mpcmulx(const mp_complex& a, const mp_complex& b,
			 mp_complex& c, int prec_words)
{
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if(debug_level >= 7)
    cerr << "MPCMULX I" << endl;

  // use 3 multiply method for complex multiplication.
  // This method requires more additions, but with this library,
  // multiplications are more expensive.
  mp_real temp1, temp2, temp4;
  mp_real::mpmulx(a.real, b.real, temp1, prec_words);
  mp_real::mpmulx(a.imag, b.imag, temp2, prec_words);
  mp_real::mpsub(temp1, temp2, temp4, prec_words);
  mp_real::mpadd(temp2, temp1, temp2, prec_words);
  mp_real temp3;
  mp_real::mpadd(a.real, a.imag, temp1, prec_words);
  mp_real::mpadd(b.real, b.imag, temp3, prec_words);
  mp_real::mpmulx(temp1, temp3, temp3, prec_words);
  mp_real::mpsub(temp3, temp2, c.imag, prec_words);
  mp_real::mpeq(temp4, c.real, prec_words);
  return;
}


void mp_complex::mpcdiv(const mp_complex& a, const mp_complex& b,
			mp_complex& c, int prec_words)
{
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if(debug_level >= 8)
    cerr << "MPCDIV I" << endl;

  if(b.real[1] == 0 && b.imag[1] == 0) {
    if(MPKER[16] != 0) {
      cerr <<"*** MPCDIV: Divisor is zero" << endl;
      error_no = 16;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  mp_real temp1, temp2, temp3, creal, denom;
  mp_real::mpmul(a.real, b.real, temp1, prec_words);
  mp_real::mpmul(a.imag, b.imag, temp2, prec_words);
  //compute denominator and do not invert.
  mp_real::mpmul(b.real, b.real, temp3, prec_words);
  mp_real::mpmul(b.imag, b.imag, denom, prec_words);
  mp_real::mpadd(temp3, denom, denom, prec_words);
  
  //real part of result
  mp_real::mpadd(temp1, temp2, temp3, prec_words);
  mp_real::mpdiv(temp3, denom, creal, prec_words);
  
  //imag part of result
  mp_real::mpsub(temp2, temp1, temp2, prec_words);
  mp_real::mpadd(a.real, a.imag, temp1, prec_words);
  mp_real::mpsub(b.real, b.imag, temp3, prec_words);
  mp_real::mpmul(temp1, temp3, temp3, prec_words);
  mp_real::mpadd(temp2, temp3, temp3, prec_words);
  mp_real::mpdiv(temp3, denom, c.imag, prec_words);
  mp_real::mpeq(creal, c.real, prec_words);
}

void mp_complex::mpcdivx(const mp_complex& a, const mp_complex& b,
			 mp_complex& c, int prec_words)
{
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(c);
    return;
  }
  if(debug_level >= 8)
    cerr << "MPCDIV I" << endl;

  if(b.real[1] == 0 && b.imag[1] == 0) {
    if(MPKER[18] != 0) {
      cerr <<"*** MPCDIVX: Divisor is zero" << endl;
      error_no = 18;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  //Check if precision level justifies advanced routine.
  int ncr = 1 <<mpmcrx;
  if(prec_words < ncr) {
    mpcdiv(a, b, c, prec_words); return;
  }

  mp_real temp1, temp2, temp3, creal, denom;
  mp_real::mpmulx(a.real, b.real, temp1, prec_words);
  mp_real::mpmulx(a.imag, b.imag, temp2, prec_words);
  //compute denominator and invert.
  mp_real::mpsqx(b.real, temp3, prec_words);
  mp_real::mpsqx(b.imag, denom, prec_words);
  mp_real::mpadd(temp3, denom, temp3, prec_words);
  mp_real f(1.0, 8);
  mp_real::mpdivx(f, temp3, denom, prec_words);
  
  //real part of result
  mp_real::mpadd(temp1, temp2, temp3, prec_words);
  mp_real::mpmulx(temp3, denom, creal, prec_words);
  
  //imag part of result
  mp_real::mpsub(temp1, temp2, temp2, prec_words);
  mp_real::mpadd(a.real, a.imag, temp1, prec_words);
  mp_real::mpsub(b.real, b.imag, temp3, prec_words);
  mp_real::mpmulx(temp1, temp3, temp3, prec_words);
  mp_real::mpsub(temp3, temp2, temp3, prec_words);
  mp_real::mpmulx(temp3, denom, c.imag, prec_words);
  mp_real::mpeq(creal, c.real, prec_words);
}

void mp_complex::mpceq(const mp_complex& a, mp_complex& b, int prec_words)
{
  mp_real::mpeq(a.real, b.real, prec_words);
  mp_real::mpeq(a.imag, b.imag, prec_words);
} 

void mp_complex::mpcsqx(const mp_complex& a, mp_complex& c, int prec_words)
{
  // use 3 multiply method for complex multiplication.
  // This method requires more additions, but with this library,
  // multiplications are more expensive.
  mp_real temp1, temp2;
  mp_real temp3;
  mp_real::mpsqx(a.real, temp1, prec_words);
  mp_real::mpsqx(a.imag, temp2, prec_words);
  mp_real::mpadd(a.real, a.imag, temp3, prec_words);
  mp_real::mpsub(temp1, temp2, c.real, prec_words);
  mp_real::mpadd(temp2, temp1, temp2, prec_words);
  mp_real::mpsqx(temp3, temp1, prec_words);
  mp_real::mpsub(temp1, temp2, c.imag, prec_words);
  return;
}

void mp_complex::mpcagx(mp_complex& a, mp_complex& b)
{
  /**
   * This performs the arithmetic-geometric mean (AGM iterations.
   * This routine is called by MPANGX.  It is not intended to be called
   * directly by the uesr.
   */
  int prec_words = mp::prec_words;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(a);
    zero(b);
    return;
  }
  int n5 = prec_words+5;
  mp_complex sk0(1.0, 1.0, n5), sk1(1.0, 1.0, n5);
  int l1 = 0;
  sk0.real[1] = 1.0;
  sk0.real[2] = 0.0;
  double s1 = 2.0;
  
  //While condition checks for convergence
  while(l1<500 && sk0.real[1] != 0.0 && 
	(sk0.real[2] < s1 || sk0.real[2] >= -2)) {
    l1++;
    s1 = sk0.real[2];
    mpcadd(a, b, sk0, prec_words);
    mpcmuld(sk0, 0.5, 0, sk1, prec_words);
    mpcmulx(a, b, sk0, prec_words);
    mpcsqrtx(sk0, b);
    mpceq(sk1, a, prec_words);
    mp_real::mpsub(a.real, b.real, sk0.real, prec_words);
  }
  if(l1 >= 500) {
    if(MPKER[12] != 0) {
      cerr <<"*** MPCAGX: Iteration limit exceeded." << endl;
      error_no = 12;
      if(MPKER[error_no] == 2) mpabrt();
    } 
    return;
  }
  if(debug_level >= 6) {
    cerr << "MPCAGX: Iter = " << l1 << endl;
    cerr << "MPCAGX: Tol. Achieved = " << sk0.real[2] << endl;
  }
  return;
}

void mp_complex::mpcsqrtx(const mp_complex& a, mp_complex& b)
{
  /**
   * This routine computes the complex square root of the 
   * MPC number C. L is the offset between real and imaginary parts
   * in A and B.  L must be at least prec_words+4.  
   */
  int prec_words = mp::prec_words;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
  }
  if(debug_level >= 5) cerr <<"MPCSQRTX I" << endl;

  if(a.real[1] == 0.0 && a.imag[1] == 0.0) {
    zero(b);
    return;
  }

  int n5 = prec_words + 5;
  mp_real sk0(0.0, n5), sk1(0.0, n5), sk2(0.0, n5);
  
  mp_real::mpsqx(a.real, sk0, prec_words);
  mp_real::mpsqx(a.imag, sk1, prec_words);
  mp_real::mpadd(sk0, sk1, sk2, prec_words);
  mp_real::mpsqrtx(sk2, sk0, prec_words);
  mp_real::mpeq(a.real, sk1, prec_words);
  sk1[1] = std::abs(sk1[1]); // take abs value.
  mp_real::mpadd(sk0, sk1, sk2, prec_words);
  mp_real::mpmuld(sk2, 0.5, 0, sk1, prec_words);
  mp_real::mpsqrtx(sk1, sk0, prec_words);
  mp_real::mpmuld(sk0, 2.0, 0, sk1, prec_words);
  if(a.real[1] >= 0.0) {
    mp_real::mpeq(sk0, b.real, prec_words);
    mp_real::mpdivx(a.imag, sk1, b.imag, prec_words);
  } else {
    mp_real::mpdivx(a.imag, sk1,  b.real, prec_words);
    b.real[1] = std::abs(b.real[1]);
    mp_real::mpeq(sk0, b.imag, prec_words);
    b.imag[1] = sign(b.imag[1], a.imag[1]);
  }

  if(debug_level >= 5) cerr <<"MPCSQRTX 0" << endl;
  return;
}

void mp_complex::mpcpwx(const mp_complex& a, int n, mp_complex& b)
{
  /**
   * This computes the N-th power of the MPC number A and places
   * the MPC result in B.  When N is zero, 1 is returned.  When N is negative,
   * the reciprocal of A ^ |N| is returned.  Before calling MPNPWX, the
   * arrays of mpuu1 and uu2 must be initialized by calling mpinix.
   * For modest levels of precision, use MPNWPR.
   *
   * This routine emplys the normal binary method for exponentiation.
   */
  const double cl2 = 1.4426950408889633;
  int prec_words = mp::prec_words;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 6) cerr << "MPCPWX I" << endl;

  int ncr = 1 << mpmcrx;
  int na1 = std::min(int(std::abs(a.real[1])), prec_words);
  int na2 = std::min(int(std::abs(a.imag[1])), prec_words);
  
  // Check if precision level of A is too low to justify the advanced routine.
  
  if(prec_words < ncr || (n > 0 && na1 < ncr && na2 < ncr)) {
    mpcpwr(a, n, b);
    return;
  }
  if(na1 == 0 && na2 == 0) {
    if(n >= 0) {
      zero(b);
      return;
    } else {
      // error : 0^N, N <= 0.
      if(MPKER[26] != 0) {
	cerr <<"*** MPCPWX: argument is zero and N is negative or zero." << endl;
	error_no = 26;
	if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }

  int n5 = prec_words+5;
  int nn = std::abs(n);
  mp_complex sk0(0.0, 0.0, n5), sk1(0.0, 0.0, n5), f1(1.0, 0.0, 8);

  if(nn == 0) {
    mpceq(f1, b, prec_words);
    return;
  } else if(nn == 1) {
    mpceq(a, b, prec_words);
  } else if(nn == 2) {
    mpcsqx(a, b, prec_words);
  } else {
    // full binary exponentiation.
    // Determine the least integer mn such that 2 & mn > nn.
    int mn = int(ANINT(cl2 * log(double(nn)) + 1.0 + mprxx));
    mpceq(a, sk0, prec_words);
    mpceq(f1, b, prec_words);
    int kn = nn, j, kk;

    for(j=1; j <= mn; j++) {
      kk = kn / 2;
      if(kn != 2 * kk) {
	mpcmulx(b, sk0, b, prec_words);
      }
      kn = kk;
      if(j < mn)
	mpcsqx(sk0, sk0, prec_words);
    }
  }
  
  // Compute the reciprocal if N is negative.
  if(n < 0) {
    mpcdivx(f1, b, sk0, prec_words);
    mpceq(sk0, b, prec_words);
  }
  return;
}

void mp_complex::mpcpwr(const mp_complex& a, int n, mp_complex& b)
{
  /**
   * This computes the N-th power of the MPC number A and places
   * the MPC result in B.  When N is zero, 1 is returned.  When N is negative,
   * the reciprocal of A ^ |N| is returned.  
   *
   * This routine employs the normal binary method for exponentiation.
   */
  int prec_words = mp::prec_words;
  
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if(debug_level >= 7) cerr << "MPCPWR I" << endl;

  int na1 = std::min(int(std::abs(a.real[1])), prec_words);
  int na2 = std::min(int(std::abs(a.imag[1])), prec_words);
  
  // Check if precision level of A is too low to justify the advanced routine.
  
  if(na1 == 0 && na2 == 0) {
    if(n >= 0) {
      zero(b);
      return;
    } else {
      // error : 0^N, N <= 0.
      if(MPKER[25] != 0) {
	cerr <<"*** MPCPWX: argument is zero and N is negative or zero." << endl;
	error_no = 25;
	if(MPKER[error_no] == 2) mpabrt();
      }
      return;
    }
  }

  int nws = prec_words;
  prec_words++;
  int n5 = prec_words+5;
  int nn = std::abs(n);
  mp_complex sk0(0.0, 0.0, n5), sk1(0.0, 0.0, n5), f1(1.0, 0.0, 8);

  if(nn == 0) {
    mpceq(f1, b, prec_words);
    prec_words = nws;
    return;
  } else if(nn == 1) {
    mpceq(a, b, prec_words);
  } else if(nn == 2) {
    mpcsqx(a, b, prec_words);
  } else {
    // full binary exponentiation.
    mpceq(f1, b, prec_words); // b = 1.0;
    mpceq(a, sk0, prec_words);
    int kn;

    for(kn = nn; kn; kn /= 2) {
      if(kn & 0x1) {
	mpcmul(b, sk0, b, prec_words);
      }
      if(kn > 1) {//dont square last iteration
	mpcsqx(sk0, sk0, prec_words);
      }
    } 
  }
  
  // Compute the reciprocal if N is negative.
  if(n < 0) {
    mpcdiv(f1, b, sk0, prec_words);
    mpceq(sk0, b, prec_words);
  }
  prec_words = nws;
  mp_real::mproun(b.real);
  mp_real::mproun(b.imag);
  return;
}

void mp_complex::mpcsqrt(const mp_complex &a, mp_complex& b)
{
  /**
   * This routine computes the complex square root of the MPC number A.
   * for extra high levels of precision, use MPCSQRTX.
   * the last word of the result is not reliable. 
   *
   *
   * This routine uses the following formula, where A1 and A2 are
   * the real and imaginary parts of A, and where R = Sqrt [A1 ^2 + A2 ^ 2]:
   *
   * B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]
   *
   * If the imaginary part of A is < 0, then the imaginary part of B
   * is also set to be < 0.
   */
  int prec_words = mp::prec_words;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    zero(b);
  }
  if(debug_level >= 6) cerr << "MPCSQRT I" << endl;
  if(a.real[1] == 0 && a.imag[1] == 0) {
    zero(b);
    return;
  }
  int n6 = prec_words+6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  
  mp_real::mpmul(a.real, a.real, sk0, prec_words);
  mp_real::mpmul(a.imag, a.imag, sk1, prec_words);
  mp_real::mpadd(sk0, sk1, sk2, prec_words);
  mp_real::mpsqrt(sk2, sk0, prec_words);
  mp_real::mpeq(a.real, sk1, prec_words);
  sk1[1] = std::abs(sk1[1]);
  mp_real::mpadd(sk0, sk1, sk2, prec_words);
  mp_real::mpmuld(sk2, 0.5, 0, sk1, prec_words);
  mp_real::mpsqrt(sk1, sk0, prec_words);
  mp_real::mpmuld(sk0, 2.0, 0, sk1, prec_words);
  if(a.real[1] >= 0) {
    mp_real::mpeq(sk0, b.real, prec_words);
    mp_real::mpdiv(a.imag, sk1, b.imag, prec_words);
  } else {
    mp_real::mpdiv(a.imag, sk1, b.real, prec_words);
    b.real[1] = std::abs(b.real[1]);
    mp_real::mpeq(sk0, b.imag, prec_words);
    b.imag[1] = sign(b.imag[1], a.imag[1]);
  }
  return;
}

mp_complex_temp exp(const mp_complex& a) {
  mp_real b;
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpexpx(a.real, mp_real::_pi, mp_real::_log2, b);
  mp_real::mpcssx(a.imag, mp_real::_pi, c.real, c.imag);
  mp_real::mpmulx(b, c.real, c.real, prec_words);
  mp_real::mpmulx(b, c.imag, c.imag, prec_words);
  return c.toTempAndDestroy();
}

mp_complex_temp log(const mp_complex& a) {
  mp_real b, c;
  mp_complex ret;
  int prec_words = mp::prec_words;
  mp_real::mpsqx(a.real, b, prec_words);
  mp_real::mpsqx(a.imag, c, prec_words);
  mp_real::mpadd(b, c, c, prec_words);
  mp_real::mplogx(c, mp_real::_pi, mp_real::_log2, b, prec_words);
  mp_real::mpmuld(b, 0.5, 0, ret.real, prec_words);
  mp_real::mpangx(a.real, a.imag, mp_real::_pi, ret.imag);
  return ret.toTempAndDestroy();
}

mp_complex_temp sin(const mp_complex &a) {
  mp_real mpt1, mpt2, mpt3, mpt4, mpt5, mpt6;
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpeq(a.imag, mpt2, prec_words);
  mpt2[1] = - mpt2[1];
  mp_real::mpexpx(mpt2,  mp_real::_pi, mp_real::_log2, mpt1);
  mpt3[1] = 1.0; mpt3[2] = 0.0; mpt3[3] = 1.0; mpt3[4] = 0.0;
  mp_real::mpdivx(mpt3, mpt1, mpt2, prec_words);
  mp_real::mpcssx(a.real, mp_real::_pi, mpt3, mpt4);
  mp_real::mpadd(mpt1, mpt2, mpt5, prec_words);
  mp_real::mpmuld(mpt5, 0.5, 0, mpt6, prec_words);
  mp_real::mpmulx(mpt6, mpt4, c.real, prec_words);
  mp_real::mpsub(mpt1, mpt2, mpt5, prec_words);
  mp_real::mpmuld(mpt5, -0.5, 0, mpt6, prec_words);
  mp_real::mpmulx(mpt6, mpt3, c.imag, prec_words);
  return c.toTempAndDestroy();
}

mp_complex_temp cos(const mp_complex &a) {
  mp_real mpt1, mpt2, mpt3, mpt4, mpt5, mpt6;
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpeq(a.imag, mpt2, prec_words);
  mpt2[1] = -mpt2[1];
  mp_real::mpexpx(mpt2, mp_real::_pi, mp_real::_log2, mpt1);
  mpt3[1] = 1.0; mpt3[2] = 0.0; mpt3[3] = 1.0; mpt3[4] = 0.0;
  mp_real::mpdivx(mpt3, mpt1, mpt2, prec_words);
  mp_real::mpcssx(a.real, mp_real::_pi, mpt3, mpt4);
  mp_real::mpadd(mpt1, mpt2, mpt5, prec_words);
  mp_real::mpmuld(mpt5, 0.5, 0, mpt6, prec_words);
  mp_real::mpmulx(mpt6, mpt3, c.real, prec_words);
  mp_real::mpsub(mpt1, mpt2, mpt5, prec_words);
  mp_real::mpmuld(mpt5, 0.5, 0, mpt6, prec_words);
  mp_real::mpmulx(mpt6, mpt4, c.imag, prec_words);
  return c.toTempAndDestroy();
}

mp_complex_temp sqr(const mp_complex &a) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_complex::mpcsqx(a, c,prec_words);
  return c.toTempAndDestroy();
}

mp_complex_temp sqrt(const mp_complex &a) {
  mp_complex c;
  mp_complex::mpcsqrtx(a, c);
  return c.toTempAndDestroy();
}

mp_real_temp abs(const mp_complex &a) {
  return sqrt(sqr(a.real) + sqr(a.imag));
}

mp_real_temp arg(const mp_complex &a) {
  return atan2(a.imag, a.real);
}

mp_complex_temp pow(const mp_complex& a, int n) {
  mp_complex c;
  mp_complex::mpcpwx(a, n, c);
  return c.toTempAndDestroy();
}

mp_complex_temp pow(const mp_complex& a, const mp_real& b) {
  mp_complex c;
  c = log(a);
  c *= b;
  c = exp(c);
  return c.toTempAndDestroy();
}

mp_complex_temp pow(const mp_complex& a, const mp_complex& b) {
  mp_complex c;
  c = log(a);
  c *= b;
  c = exp(c);
  return c.toTempAndDestroy();
}

