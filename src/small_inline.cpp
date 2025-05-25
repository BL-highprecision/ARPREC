#include <arprec/mp_real.h>

#ifdef ARPREC_VACPP_BUILTINS_H
/* for IBM VisualAge C++ __fmadd */
#include <builtins.h>
#endif

#define _SPLITTER 134217729.0               /* = 2^27 + 1 */

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline double quick_two_sum(double a, double b, double &err) {
  double s = a + b;
  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline double two_sum(double a, double b, double &err) {
  double s = a + b;
  double bb = s - a;
  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

#ifndef ARPREC_FMS
/* Computes high word and lo word of a */
inline void split(double a, double &hi, double &lo) {
  double temp;
  temp = _SPLITTER * a;
  hi = temp - (temp - a);
  lo = a - hi;
}
#endif

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline double two_sqr(double a, double &err) {
#ifdef ARPREC_FMS
  double p = a * a;
  err = ARPREC_FMS(a, a, p);
  return p;
#else
  double hi, lo;
  double q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Note :  The below functions use the IEEE trick 
 *   of first adding TWO_TO_THE_52, and then subtracting
 *   it again, produces a rounded (sometimes incorrectly rounded) integer
 *   value when the input is between 0 and 2^52-1. 
 *
 *  routines prefaced with SLOPPY are slightly faster,
 *  but can round strangly (may be off by one).
 */

#define TWO_TO_THE_52 4503599627370496.0
inline double SLOPPY_ANINT_POSITIVE(double a) { 
  // performs anint, with sometimes incorrect rounding.
  // Assumes that the input is in [0, 2^52).
  return (a + TWO_TO_THE_52) - TWO_TO_THE_52;
}
inline double FLOOR_POSITIVE(double a) {
  // performs floor, correctly.
  // Assumes that the input is in [0, 2^52).
  double b = (a + TWO_TO_THE_52) - TWO_TO_THE_52;
  if (b>a) return b-1.0; else return b;
}  
inline double CEIL_POSITIVE(double a) {
  // performs floor, correctly.
  // Assumes that the input is in [0, 2^52).
  double b = (a + TWO_TO_THE_52) - TWO_TO_THE_52;
  if (b<a) return b+1.0; else return b;
}  
inline double AINT(double a) {
  // performs aint, correctly.
  // Assumes that the input is in (-2^52, 2^52).
  double b = a;
  if (a>=0) {
    b = (b + TWO_TO_THE_52) - TWO_TO_THE_52;
    if(b>a) return b-1.0; else return b;
  }
  else {
    b = (b - TWO_TO_THE_52) + TWO_TO_THE_52;
    if(b<a) return b+1.0; else return b;
  }
}  
inline double POSITIVE_AINT(double a) {
  // performs aint, correctly.
  // Assumes that the input is in [0, 2^52).
  double b = a;
  b = (b + TWO_TO_THE_52) - TWO_TO_THE_52;
  if (b > a) return b - 1.0; else return b;
}
inline double ANINT(double a) { 
  // performs anint, correctly.
  // Assumes that the input is in (-2^52, 2^52).
  return a >= 0.0 ? AINT(a + 0.5) : AINT(a - 0.5); 
}
inline double SLOPPY_ANINT(double a) {
  // this one is correct most of the time.
  // performs anint, with possible rounding errors.
  // Assumes that the input is in (-2^52, 2^52).
  a = (a + TWO_TO_THE_52) - TWO_TO_THE_52; //for positives.
  a = (a - TWO_TO_THE_52) + TWO_TO_THE_52; //for negatives.
  return a;
}
#undef TWO_TO_THE_52

inline void dd_add_d(double a[], double b, double c[]) {
  // This routine adds double-double A with double B, where B
  // corresponds to the low word of A.
  //
  // The algorithm is modified from Briggs' DD + DD.
  //
  double s1, s2, t1, t2;

  t1 = two_sum(a[1], b, t2);
  s1 = quick_two_sum(a[0], t1, s2);
  s2 += t2;
  t1 = quick_two_sum(s1, s2, t2);
  c[0] = t1;
  c[1] = t2;
}

inline void dd_add_dd(double a[], double b[], double c[]) {
  // This routine adds double-double A with double-double B.
  //
  // The algorithm due to D. Bailey.
  //
  double s, e;

  s = two_sum(a[0], b[0], e);
  e += a[1];
  e += b[1];
  s = quick_two_sum(s, e, e);
  c[0] = s;
  c[1] = e;
}

inline double mp_two_prod(double a, double b, double &err) {
  // This demonstrates the scheme to multiply two DP ARPREC entries, each
  // representing a 52-bit mantissa integer, and produce the double-DP result.
  double s1;
#ifdef ARPREC_FMS
  double p = a * b;
  err = ARPREC_FMS(a, b, p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  double p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  s1 = p;
#endif
  //assert(a < pow(2.0, 52) && a > -pow(2.0, 52));
  //assert(b < pow(2.0, 52) && b > -pow(2.0, 52));

  // This clean-up makes sure that each word has at most mpnbt bits
  // The final product is 2^52 * s1 + err
  double t = s1 * mp_real::mprdx;
  s1 = SLOPPY_ANINT (t);
  err += mp_real::mpbdx * (t - s1);

  return s1;
}

inline double mp_two_prod_positive(double a, double b, double &err) {
  // This demonstrates the scheme to multiply two DP ARPREC entries, each
  // representing a 52-bit mantissa integer, and produce the double-DP result.
  // This assumes that a, b >=0.0
  // This routine is intended for use with mpmul and mpmuld only.
  double s1;
#ifdef ARPREC_FMS
  double p = a * b;
  err = ARPREC_FMS(a, b, p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  double p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  s1 = p;
#endif
  //assert(a < pow(2.0, 52) && a > -pow(2.0, 52));
  //assert(b < pow(2.0, 52) && b > -pow(2.0, 52));

  // This clean-up makes sure that each word has at most mpnbt bits
  // The final product is 2^52 * s1 + err
  double t = s1 * mp_real::mprdx;
  s1 = SLOPPY_ANINT_POSITIVE(t); // ok since inputs always positive
  err += mp_real::mpbdx * (t - s1);

  return s1;
}
