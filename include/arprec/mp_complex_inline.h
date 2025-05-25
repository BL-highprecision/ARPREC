/*
 * include/arprec/mp_complex_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef MP_COMPLEX_INLINE_H
#define MP_COMPLEX_INLINE_H

inline mp_complex_temp operator+(const mp_complex &a, const mp_complex &b) {
  return mp_complex_temp(a.real+b.real, a.imag+b.imag);
}

inline  mp_complex_temp operator+(const mp_complex &a, const mp_real &breal) {
  mp_real temp(a.imag);
  return mp_complex_temp(a.real+breal, temp.toTempAndDestroy());
}

inline  mp_complex_temp operator+(const mp_real &breal, const mp_complex &a) {
  mp_real temp(a.imag);
  return mp_complex_temp(a.real+breal, temp.toTempAndDestroy());
}

inline mp_complex_temp operator-(const mp_complex &a, const mp_complex &b) {
  return mp_complex_temp(a.real-b.real, a.imag-b.imag);
}

inline  mp_complex_temp operator-(const mp_complex &a, const mp_real &breal) {
  mp_real temp(a.imag);
  return mp_complex_temp(a.real - breal, temp.toTempAndDestroy());
}

inline  mp_complex_temp operator-(const mp_real &breal, const mp_complex &a) {
  mp_real temp(a.imag);
  // bug fix -- XSL
  temp[1] = -temp[1];
  return mp_complex_temp(breal - a.real, temp.toTempAndDestroy());
}

inline mp_complex_temp operator*(const mp_complex &a, const mp_complex &b) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_complex::mpcmulx(a, b, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_complex_temp operator/(const mp_complex &a, const mp_complex &b) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_complex::mpcdivx(a, b, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_complex_temp operator/(const mp_complex &a, const mp_real &b) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpdivx(a.real, b, c.real, prec_words);
  mp_real::mpdivx(a.imag, b, c.imag, prec_words);
  return c.toTempAndDestroy();
}

inline mp_complex_temp operator/(const mp_complex &a, double b) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpdivd(a.real, b, 0, c.real, prec_words);
  mp_real::mpdivd(a.imag, b, 0, c.imag, prec_words);
  return c.toTempAndDestroy();
}

inline mp_complex_temp operator/(const mp_real &a, const mp_complex &b) {
  mp_real c;
  c = sqr(b.real) + sqr(b.imag);
  c = a/c;
  mp_real_temp ret_imag(b.imag * c);
  ret_imag.mpr[1] = -ret_imag.mpr[1]; //conjugate result.
  return mp_complex_temp(b.real*c,ret_imag); 
}

inline mp_complex_temp operator/(double b, const mp_complex& a) {
  return mp_real(b) / a;
}

inline mp_complex_temp operator/(int b, const mp_complex& a) {
  return mp_real(static_cast<double>(b)) / a;
}

inline mp_complex_temp operator*(const mp_complex& a, double b) {
  mp_complex ret;
  int prec_words = mp::prec_words;
  mp_complex::mpcmuld(a, b, 0, ret, prec_words);
  return ret.toTempAndDestroy();
}

inline mp_complex_temp operator*(double b, const mp_complex& a) {
  return a * b;
}

inline mp_complex_temp operator*(const mp_complex& a, const mp_real& b) {
  mp_complex c;
  int prec_words = mp::prec_words;
  mp_real::mpmulx(a.real, b, c.real, prec_words);
  mp_real::mpmulx(a.imag, b, c.imag, prec_words);
  return c.toTempAndDestroy();
}

inline mp_complex_temp operator*(const mp_real& b, const mp_complex& a) {
  return a * b;
}

inline mp_complex &mp_complex::operator+=(const mp_complex &other) {
  int prec_words = mp::prec_words;
  mp_real::mpadd(this->real, other.real, this->real, prec_words);
  mp_real::mpadd(this->imag, other.imag, this->imag, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator-=(const mp_complex &other) {
  int prec_words = mp::prec_words;
  mp_real::mpsub(this->real, other.real, this->real, prec_words);
  mp_real::mpsub(this->imag, other.imag, this->imag, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator*=(const mp_complex &other) {
  mp_complex temp;
  int prec_words = mp::prec_words;
  mpcmulx(*this, other, temp, prec_words);
  mpceq(temp, *this, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator/=(const mp_complex &other) {
  mp_complex temp;
  int prec_words = mp::prec_words;
  mpcdivx(*this, other, temp, prec_words);
  mpceq(temp, *this, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator+=(const mp_real &other) {
  int prec_words = mp::prec_words;
  mp_real::mpadd(this->real, other, this->real, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator-=(const mp_real &other) {
  int prec_words = mp::prec_words;
  mp_real::mpsub(this->real, other, this->real, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator*=(const mp_real &other) {
  int prec_words = mp::prec_words;
  mp_real::mpmulx(this->real, other, this->real, prec_words);
  mp_real::mpmulx(this->imag, other, this->imag, prec_words);
  return *this;
}

inline mp_complex &mp_complex::operator/=(const mp_real &other) {
  int prec_words = mp::prec_words;
  mp_real::mpdivx(this->real, other, this->real, prec_words);
  mp_real::mpdivx(this->imag, other, this->imag, prec_words);
  return *this;
}

inline bool operator==(const mp_complex &a, const mp_complex &b) {
  return (a.real == b.real) && (a.imag == b.imag);
}

inline bool operator!=(const mp_complex &a, const mp_complex &b) {
  return (a.real != b.real) || (a.imag != b.imag);
}

inline mp_complex &mp_complex::operator=(const mp_complex &other) {
  this->real = other.real;
  this->imag = other.imag;
  return *this;
}

inline mp_complex &mp_complex::operator=(const mp_real &other) {
  this->real = other;
  mp_real::zero(this->imag);
  return *this;
}

inline mp_complex &mp_complex::operator=(mp_complex_temp &other) {
  this->real = other.real;
  this->imag = other.imag;
  return *this;
}

inline mp_complex &mp_complex::operator/=(double other) {
  return *this = *this / other;
}

inline mp_complex &mp_complex::operator/=(int other) {
  return *this = *this / other;
}

#endif /* MP_COMPLEX_INLINE_H */
