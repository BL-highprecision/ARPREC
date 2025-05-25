/*
 * include/arprec/mp_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MP_INLINE_H
#define __MP_INLINE_H

#include <cstring>
#include <arprec/mp.h>
#include <arprec/double_inline.h>

/* mp-real + mp-real */
inline mp_real_temp operator+(const mp_real &a, const mp_real &b) {
  mp_real c;
  mp_real::mpadd(a, b, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* mp-real - mp-real */
inline mp_real_temp operator-(const mp_real &a, const mp_real &b) {
  mp_real c;
  mp_real::mpsub(a, b, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* mp-real * mp-real */
inline mp_real_temp operator*(const mp_real &a, const mp_real &b) {
  mp_real c;
  mp_real::mpmulx(a, b, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* mp-real * double */
inline mp_real_temp operator*(const mp_real &a, const double b) {
  mp_real c;
  mp_real::mpmuld(a, b, 0, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* double * mp-real */
inline mp_real_temp operator*(const double b, const mp_real &a) {
  mp_real c;
  mp_real::mpmuld(a, b, 0, c, mp::prec_words);
  return c.toTempAndDestroy();
}


/* mp-real / mp-real */
inline mp_real_temp operator/(const mp_real &a, const mp_real &b) {
  mp_real c;
  mp_real::mpdivx(a, b, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* mp-real / double */
inline mp_real_temp operator/(const mp_real &a, const double b) {
  mp_real c;
  mp_real::mpdivd(a, b, 0, c, mp::prec_words);
  return c.toTempAndDestroy();
}

/* double / mp_real */
inline mp_real_temp operator/(const double b, const mp_real &a) {
  return mp_real(b, 8) / a;
}

inline mp_real& mp_real::operator=(const int& ib) {
  if(ib == 0) { zero(*this); return *this; }
  mpdmc(ib, 0, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator=(const double& db) {
  int prec_words = mp::prec_words;
  if(db == 0.0) { zero(*this); return *this; }
  mpdmc(db, 0, *this, prec_words);
  return *this;
}

inline mp_real& mp_real::operator=(const mp_real& jb) {
  if(&jb != this)
    mpeq(jb, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator=(const char *a) {
  mp_real::mpdexc(a, static_cast<unsigned int>(strlen(a)), *this);
  return *this;
}

inline mp_real& mp_real::operator+=(const mp_real &b) {
  //same variable input-output safe
  mp_real::mpadd(*this, b, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator*=(const mp_real &b) {
  //same variable input-output safe
  mp_real::mpmulx(*this, b, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator*=(double b) {
  //same variable input-output safe
  mp_real::mpmuld(*this, b, 0, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator-=(const mp_real &b) {
  //same variable input-output safe (for the first & third arguments only!)
  mp_real::mpsub(*this, b, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator/=(const mp_real &b) {
  //same variable input-output safe
  mp_real::mpdivx(*this, b, *this, mp::prec_words);
  return *this;
}

inline mp_real& mp_real::operator/=(double b) {
  //same variable input-output safe
  mp_real::mpdivd(*this, b, 0, *this, mp::prec_words);
  return *this;
}

inline bool operator>(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) == 1;
}

inline bool operator>=(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) >= 0;
}

inline bool operator<(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) < 0;
}

inline bool operator<=(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) <= 0;
}

inline bool operator==(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) == 0;
}

inline bool operator!=(const mp_real& a, const mp_real& b) {
  return mp_real::mpcpr(a, b) != 0;
}

/* (- mp-real) */
inline mp_real_temp operator-(const mp_real &a) {
  mp_real b(a);
  b[1] = -b[1];
  return b.toTempAndDestroy();
}

inline mp_real_temp abs(const mp_real& a) {
  mp_real ret;
  mp_real::mpeq(a, ret, mp::prec_words);
  ret[1] = fabs(ret[1]);
  return ret.toTempAndDestroy();
}

/**
 * Operations with int need to have explicit operators.
 * If they did not, the mp_real(int) would be used, 
 * which creates a mp_real of a certain size, rather
 * than an mp_real representing an integer.
 */
inline mp_real_temp operator+(const mp_real &a, int b) {
  return a + mp_real(double(b), 8);
}

inline mp_real_temp operator+(int b, const mp_real &a) {
  return a + mp_real(double(b), 8);
}

inline mp_real operator+=(mp_real &a, int b) {
  return a += mp_real(double(b), 8);
}

inline mp_real_temp operator-(int b, const mp_real &a) {
  return mp_real(double(b), 8) - a;
}

inline mp_real_temp operator-(const mp_real &a, int b) {
  return a - mp_real(double(b), 8);
}

inline mp_real operator-=(mp_real &a, int b) {
  return a -= mp_real(double(b), 8);
}

inline mp_real_temp operator*(const mp_real &a, int b) {
  return a * double(b);
}

inline mp_real operator*=(mp_real &a, int b) {
  return a *= double(b);
}

inline mp_real_temp operator*(int b, const mp_real &a) {
  return a * double(b);
}

inline mp_real_temp operator/(const mp_real &a, int b) {
  return a / double(b);
}

inline mp_real operator/=(mp_real &a, int b) {
  return a /= double(b);
}

inline mp_real_temp operator/(int b, const mp_real &a) {
  return mp_real(b, 8) / a;
}

inline bool operator>(int b, const mp_real &a) {
  return mp_real(double(b)) > a;
}

inline bool operator>(const mp_real &a, int b) {
  return a > mp_real(double(b));
}

inline bool operator<(int b, const mp_real &a) {
  return mp_real(double(b)) < a;
}

inline bool operator<(const mp_real &a, int b) {
  return a < mp_real(double(b));
}

inline bool operator<=(int b, const mp_real &a) {
  return mp_real(double(b)) <= a;
}

inline bool operator<=(const mp_real &a, int b) {
  return a <= mp_real(double(b));
}

inline bool operator>=(int b, const mp_real &a) {
  return mp_real(double(b)) >= a;
}

inline bool operator>=(const mp_real &a, int b) {
  return a >= mp_real(double(b));
}

inline bool operator==(int b, const mp_real &a) {
  return mp_real(double(b)) == a;
}

inline bool operator==(const mp_real &a, int b) {
  return a == mp_real(double(b));
}

inline bool operator!=(int b, const mp_real &a) {
  return mp_real(double(b)) != a;
}

inline bool operator!=(const mp_real &a, int b) {
  return a != mp_real(double(b));
}

/**
 * operators with double
 */

inline mp_real_temp operator+(const mp_real &a, double b) {
  return a + mp_real(b, 8);
}

inline mp_real_temp operator+(double b, const mp_real &a) {
  return a + mp_real(b, 8);
}

inline mp_real operator+=(mp_real &a, double b) {
  return a += mp_real(double(b), 8);
}

inline mp_real_temp operator-(double b, const mp_real &a) {
  return mp_real(b, 8) - a;
}

inline mp_real_temp operator-(const mp_real &a, double b) {
  return a - mp_real(b, 8);
}

inline mp_real operator-=(mp_real &a, double b) {
  return a -= mp_real(b, 8);
}

inline bool operator>(double b, const mp_real &a) {
  return mp_real(b) > a;
}

inline bool operator>(const mp_real &a, double b) {
  return a > mp_real(b);
}

inline bool operator<(double b, const mp_real &a) {
  return mp_real(b) < a;
}

inline bool operator<(const mp_real &a, double b) {
  return a < mp_real(b);
}

inline bool operator<=(double b, const mp_real &a) {
  return mp_real(b) <= a;
}

inline bool operator<=(const mp_real &a, double b) {
  return a <= mp_real(b);
}

inline bool operator>=(double b, const mp_real &a) {
  return mp_real(b) >= a;
}

inline bool operator>=(const mp_real &a, double b) {
  return a >= mp_real(b);
}

inline bool operator==(double b, const mp_real &a) {
  return mp_real(b) == a;
}

inline bool operator==(const mp_real &a, double b) {
  return a == mp_real(b);
}

inline bool operator!=(double b, const mp_real &a) {
  return mp_real(b) != a;
}

inline bool operator!=(const mp_real &a, double b) {
  return a != mp_real(b);
}

inline void mp_real::reallocate(unsigned int sz) {
  if (mpr && alloc) 
    delete [] mpr;
  allocate(sz);
}

inline mp_real::mp_real() {
  if (prec_words > 0) 
    allocate(static_cast<unsigned int>(prec_words + 5));
  else 
    allocate(0U);
}

inline mp_real::mp_real(const mp_real &x) {
  unsigned int sz = static_cast<unsigned int>(x[0]);
  allocate(sz);
  std::memcpy(mpr, x.mpr, sizeof(double) * sz);
}

inline mp_real::mp_real(double a, int sz) {
  allocate(sz);
  if (sz) mpdmc(a, 0, *this, prec_words);
}

inline mp_real::mp_real(const mp_real_temp& x) {
  alloc = true;
  mpr = x.mpr;
}

inline mp_real::mp_real(const double *a) {
  alloc = false;
  mpr = const_cast<double *>(a);
}

inline mp_real::mp_real(double *a) {
  alloc = false;
  mpr = a;
}

#endif
