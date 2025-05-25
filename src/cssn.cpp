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

using std::cout;

void mp_real::mpcssn(const mp_real& a, const mp_real& pi, mp_real &x,
            mp_real &y, int prec_words)
{
  /**
   * This computes the cosine and sine of the MP number A and returns
   * the two MP results in X and Y, respectively.  Pi is the MP value
   * of Pi computed by a previous call to MPPI.  For extra high
   * levels of precision, use MPCSSX.  The last word of the result
   * is not reliable.  Debug output starts with debug_level == 6.
   *
   * This routine uses the conventional Taylor's series for Sin (s) :
   *  
   * Sin (s) = s - s^3 / 3! + s ^5 / 5! - s^7 / 7! ...
   *  
   * where s = t - a * pi / 2 - b * pi / 16 and the integers a and b
   * are chosen to minimize the absolute value of s.  We can then
   * compute:
   * 
   * Sin (t) = Sin (s + a * pi / 2 + b * pi / 256)
   * Cos (t) = Cos (s + a * pi / 2 + b * pi / 256)
   * 
   * by applying elementary trig identities for sums.  The sine and
   * cosine of b * pi / 16 are of the form
   *
   *        1/2 * Sqrt {2 +- Sqrt [2 +- Sqrt(2)]}.
   *
   * Reducing t in this manner insures that -Pi / 512 < s <= Pi / 512,
   * which accelerates convergence in the above series.
   *
   * The sines and cosines for values (b * pi / 256) where b in
   * an integer are stored in the two tables below.
   */
  static mp_real *pi_over_256_sine_table[129] = 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0};
  static mp_real *pi_over_256_cosine_table[129] = 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0};
  const double cpi = 3.141592653589793;
  int nq = 4;
  double t1, t2;
  int n1, na, neg=0;
  
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(x); zero(y); return;
  }
  if (debug_level >= 6) cout << "\nMPVSSN I";

  na = std::min(int(std::abs(a[1])), prec_words);
  if(na == 0) {
    //x = 1.0, y = zero.
    x[1] = 1.0; x[2] = 0.0; x[3] = 1.0;
    zero(y);
    return;
  }

  // Check is Pi has been precomputed.
  
  mpmdc (pi, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - cpi) > mprx2) {
    if(MPKER[28] != 0) {
      cout <<"\n***MPCSSN: PI must be precomputed.";
      error_no = 28;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }

  int nws = prec_words;
  prec_words++;
  int n5 = prec_words+5;
  int ka, kb, kc;
  mp_real sk0(0.0, n5), sk1(0.0, n5), sk2(0.0, n5), sk3(0.0, n5);
  mp_real sk4(0.0, n5), sk5(0.0, n5), sk6(0.0, n5), f(1.0, 6);

  //        Reduce to between -Pi and Pi.

  mpmuld(pi, 2.0, 0, sk0, prec_words);
  mpdiv(a, sk0, sk1, prec_words);
  mpnint(sk1, sk2, prec_words);
  mpsub(sk1, sk2, sk3, prec_words);
  
  // Determine rearest multiple of Pi / 2, and within a quadrant, the
  // nearest multiple of Pi / 256.  Through most of the rest of this
  // subroutine, KA and KB are the integers a and b of the algorithm
  // above.

  mpmdc(sk3, t1, n1, prec_words);
  if(n1 >= -mpnbt) {
    t1 = ldexp(t1, n1);
    t2 = 4.0 * t1;
    t1 = anint (double(t2));
    ka = int(t1);
    kb = int(anint (128.0 * (t2 - ka)));
  } else {
    ka = 0;
    kb = 0;
  }
  t1 = (128 * ka + kb) / 512.0;
  mpdmc(t1, 0, sk1, prec_words);
  mpsub(sk3, sk1, sk2, prec_words);
  mpmul(sk2, sk0, sk1, prec_words);

  // compute consine and sine of the reduced argument s.

  if(sk1[1] == 0.0) { // if sk1 == zero
    zero(sk0);
    nq = 0;
  } else {
    //Divide by 2^nq (possibly), fix after series has converged.
    if(sk1.mpr[2] < -1 || sk1.mpr[3] < mpbdx/4096.0) {
      nq = 0;
    } else {
      mpdivd(sk1, 1.0, nq, sk1, prec_words);
    }

    // Compute the Taylor's series now.
    mpeq(sk1, sk0, prec_words);
    mpmul(sk0, sk0, sk2, prec_words);
    
    int l1=0; //iteration count.
    int term_prec;
    
    neg = sk1.mpr[1] < 0.0 ? 1 : 0;
    do {
      l1++;
      t2 = - (2.0 * l1) * (2.0 * l1 + 1.0);

        // compute this term with term_prec words of precision only.
      term_prec = std::min(nws+1, nws+int(sk2.mpr[2]+sk1.mpr[2]-sk0.mpr[2])+2);
      prec_words = std::max(0, term_prec); 
      mpmul(sk1, sk2, sk3, prec_words);
      mpdivd(sk3, t2, 0, sk1, prec_words);
      prec_words = nws+1; // full precision to add term in.
      mpadd(sk1, sk0, sk0, prec_words);
      //the above line needs to change if mpadd is not safe for 
      // same variable input/output.
      
      // Check for convergence of the series in the loop condition
    } while(l1 < 10000 &&
            (sk1[1] != 0.0 && sk1[2] >= sk0[2] - prec_words));

    if(l1 >= 10000) {
      if(MPKER[29] != 0) {
        cout <<"\n*** MPCSSN: Iteration limit exceeded.";
        error_no = 29;
        if(MPKER[error_no] == 2) mpabrt();
        prec_words = nws;
        return;
      }
    }
    //answer needs to end up in sk0.
    if(nq) {
      // Perform double angle formulas, 
      // Cos (s) = 1 - 2 * Sin^2(s/2) = 2 * Cos^2(s/2) - 1 
      mpmul(sk0, sk0, sk1, prec_words);
      mpmuld(sk1, 2.0, 0, sk2, prec_words);
      mpsub(f, sk2, sk0, prec_words);
      for(int i=1;i<nq;i++) {
        mpmul(sk0, sk0, sk1, prec_words);
        mpmuld(sk1, 2.0, 0, sk2, prec_words);
        mpsub(sk2, f, sk0, prec_words);
      }
    }      
  }

  if(nq) {
    // sk0 currently holds Cos(s).
    // Compute Sin (s) = Sqrt( 1 - Cos^2(s));
    //mpeq(sk0, sk1, prec_words);
    mpmul(sk0, sk0, sk2, prec_words);
    mpsub(f, sk2, sk3, prec_words);
    mpsqrt(sk3, sk1, prec_words);
    if(neg) {
      sk1.mpr[1] = -sk1.mpr[1];
    }
  } else {
    // sk0 currently holds Sin(s).
    // Compute Cos (s) = Sqrt( 1 - Sin^2(s));
    mpeq(sk0, sk1, prec_words);
    mpmul(sk0, sk0, sk2, prec_words);
    mpsub(f, sk2, sk3, prec_words);
    mpsqrt(sk3, sk0, prec_words);
  }    

  // Now sk0 holds Cos(s), sk1 holds Sin(s).

  // Compute cosine and sine of b * Pi / 512; or, 
  //   get it from the table.  
  
  kc = std::abs(kb);
  if(pi_over_256_sine_table[kc] && 
     (pi_over_256_sine_table[kc]->mpr[0] >= sk0.mpr[0])) {
    mpeq(*pi_over_256_cosine_table[kc], sk2, prec_words);
    mpeq(*pi_over_256_sine_table[kc], sk3, prec_words);
  } else {
    f[FST_M] = 2.0; // f = 2.0

    if(kc == 0) {
      sk2[1] = 1.0; sk2[2] = 0.0; sk2[3] = 1.0;// sk2 = 1.0;
      zero(sk3);
    } else {
      switch(kc % 8) {
      case 0:
        //sk4 = 2.0 == 2*cos(0) 
        sk4[1] = 1.0; sk4[2] = 0.0; sk4[3] = 2.0;
        break;
      case 7:
      case 1:
        mpsqrt(f, sk4, prec_words);
        mpadd(f, sk4, sk5, prec_words);
        mpsqrt(sk5, sk4, prec_words);
        break;
      case 6:
      case 2:
        mpsqrt(f, sk4, prec_words);
        break;
      case 5: 
      case 3:
        mpsqrt(f, sk4, prec_words);
        mpsub(f, sk4, sk5, prec_words);
        mpsqrt(sk5, sk4, prec_words);
        break;
      case 4:
        zero(sk4);
      }
      // if kc * Pi/8 is on the negative half of the unit circle...
      if(((kc+4)/8) & 0x1) sk4[1] = -sk4[1];
      // now sk4 holds 2 * Cos (kc * Pi / 8)
      mpadd(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk4, prec_words);
      if(((kc+8)/16) & 0x1) sk4[1] = -sk4[1];
      // now sk4 holds 2 * Cos (kc * Pi / 16)
      mpadd(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk4, prec_words);
      if(((kc+16)/32) & 0x1) sk4[1] = -sk4[1];
      // now sk4 holds 2 * Cos (kc * Pi / 32)
      mpadd(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk4, prec_words);
      if(((kc+32)/64) & 0x1) sk4[1] = -sk4[1];
      // now sk4 holds 2 * Cos (kc * Pi / 64)

      mpadd(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk4, prec_words);
      // now sk4 holds 2 * Cos (kc * Pi / 128)

      // do for all kc != 0 
      mpadd(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk3, prec_words);
      mpmuld(sk3, 0.5, 0, sk2, prec_words);
      mpsub(f, sk4, sk5, prec_words);
      mpsqrt(sk5, sk4, prec_words);
      mpmuld(sk4, 0.5, 0, sk3, prec_words);
    }
    mp_real *new_sine, *new_cosine;
    new_cosine = new mp_real(sk2);
    new_sine = new mp_real(sk3);
    if(pi_over_256_sine_table[kc]) {
      //Required precision may have increased,
      // Throw away old table members.
      delete pi_over_256_sine_table[kc];
      delete pi_over_256_cosine_table[kc];
    }
    pi_over_256_sine_table[kc] = new_sine;
    pi_over_256_cosine_table[kc] = new_cosine;
  }
  if (kb < 0) sk3[1] = -sk3[1];
  // Now sk2 holds Cos (b * Pi / 256), 
  // sk3 holds Cos (b * Pi / 256).


  // Apply the trig summation identities to compute cosine and sine
  // of s + b * Pi / 256;  

  mpmul(sk0, sk2, sk4, prec_words);
  mpmul(sk1, sk3, sk5, prec_words);
  mpsub(sk4, sk5, sk6, prec_words);
  mpmul(sk1, sk2, sk4, prec_words);
  mpmul(sk0, sk3, sk5, prec_words);
  mpadd(sk4, sk5, sk1, prec_words);
  mpeq(sk6, sk0, prec_words);
  

  // This code in effect applies the trig summation identities for
  // (s + b * Pi / 256) + a * Pi / 2.
  
  switch(ka) {
  case 0: 
    mpeq(sk0, x, prec_words);
    mpeq(sk1, y, prec_words);
    break;
  case 1:
    mpeq(sk1, x, prec_words);
    x[1] = - x[1];
    mpeq(sk0, y, prec_words);
    break;
  case -1:
    mpeq(sk1, x, prec_words);
    mpeq(sk0, y, prec_words);
    y[1] = -y[1];
    break;
  case 2:
  case -2:
    mpeq(sk0, x, prec_words);
    x[1] = -x[1];
    mpeq(sk1, y, prec_words);
    y[1] = -y[1];
    break;
  }
  
  // Restore original precision level.
  
  prec_words = nws;
  mproun(x);
  mproun(y);
  
  if(debug_level >= 6) cout << "\nMPCSSN done : sin = "<<x<<"\t cos = "<<y;

  return;
}

