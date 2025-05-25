/*
 * include/arprec/mp_int_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MPINT_INLINE_H
#define __MPINT_INLINE_H

inline bool operator==(const mp_int& a, const mp_real& b) {
  return operator==(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator!=(const mp_int& a, const mp_real& b) {
  return operator!=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>(const mp_int& a, const mp_real& b) {
  return operator>(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<(const mp_int& a, const mp_real& b) {
  return operator<(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>=(const mp_int& a, const mp_real& b) {
  return operator>=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<=(const mp_int& a, const mp_real& b) {
  return operator<=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator==(const mp_real& a, const mp_int& b) {
  return operator==(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator!=(const mp_real& a, const mp_int& b) {
  return operator!=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>(const mp_real& a, const mp_int& b) {
  return operator>(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<(const mp_real& a, const mp_int& b) {
  return operator<(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>=(const mp_real& a, const mp_int& b) {
  return operator>=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<=(const mp_real& a, const mp_int& b) {
  return operator<=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator==(const mp_int& a, const mp_int& b) {
  return operator==(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator!=(const mp_int& a, const mp_int& b) {
  return operator!=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>(const mp_int& a, const mp_int& b) {
  return operator>(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<(const mp_int& a, const mp_int& b) {
  return operator<(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator>=(const mp_int& a, const mp_int& b) {
  return operator>=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline bool operator<=(const mp_int& a, const mp_int& b) {
  return operator<=(static_cast<mp_real>(a), static_cast<mp_real>(b));
}

inline mp_int_temp operator+(const mp_int& a, const mp_int& b) {
  mp_int c;
  int prec_words = mp::prec_words;
  mp_real::mpadd(a, b, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator+(const mp_int& a, int b) {
  mp_int c, temp(b);
  int prec_words = mp::prec_words;
  mp_real::mpadd(a, temp, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator+(int b, const mp_int& a) {
  mp_int c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpadd(a, t, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_real_temp operator+(const mp_int& a, double b) {
  mp_real c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpadd(a, t, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_real_temp operator+(double b, const mp_int& a) {
  mp_real c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpadd(a, t, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator-(const mp_int& a, const mp_int& b) {
  mp_int c;
  int prec_words = mp::prec_words;
  mp_real::mpsub(a, b, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator-(const mp_int& a, int b) {
  mp_int c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpsub(a, t, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator-(int b, const mp_int& a) {
  mp_int c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpsub(t, a, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_real_temp operator-(const mp_int& a, double b) {
  mp_real c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpsub(a, t, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_real_temp operator-(double b, const mp_int& a) {
  mp_real c, t(b); 
  int prec_words = mp::prec_words;
  mp_real::mpsub(t, a, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator*(const mp_int& a, const mp_int& b) {
  mp_int c;
  int prec_words = mp::prec_words;
  mp_real::mpmulx(a, b, c, prec_words);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator*(const mp_int& a, int b) {
  mp_int c;
  int prec_words = mp::prec_words;
  mp_real::mpmuld(a, double(b), 0, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_int_temp operator*(int b, const mp_int& a) {
  mp_int c;
  int prec_words = mp::prec_words;
  mp_real::mpmuld(a, double(b), 0, c, prec_words);
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

inline mp_int& mp_int::operator/=(const mp_int& other) {
  *this = *this / other;
  mp_int::ovcheck(*this);
  return *this;
}

inline mp_int& mp_int::operator/=(int other) {
  *this = *this / other;
  mp_int::ovcheck(*this);
  return *this;
}

inline mp_int_temp operator%(const mp_int& a, const mp_int &b) {
  mp_int c(static_cast<int>(a[2] - b[2] + 11.0));
  mp_int d(static_cast<int>(a[0]));
  divrem(a, b, c, d);
  return d.toTempAndDestroy();
}

inline int operator%(int a, const mp_int &b) {
  if(b[1] == 0.0) {
    std::cerr << "\n*** MPINT, operator%(int, mp_int) : modulus by zero";
    mp_int::mpabrt();
  }
  if(b[2] > 0.0 || b[FST_M] > std::abs(a)) {
    return a;
  }
  return a % static_cast<int>(b[FST_M]);
}

inline int operator%(const mp_int& a, int b) {
  mp_int t;
  t = a / b;
  t *= b;
  t = a - t;
  //t should now have at most 1 mantissa word.
  return int(t[FST_M] * t[1]);  
}

inline mp_int& mp_int::operator%=(const mp_int& other) {
  mp_int c(static_cast<int>(mpr[2] - other[2] + 11.0));
  divrem(*this, other, c, *this);
  return *this;
}

inline mp_int& mp_int::operator%=(int other) {
  return *this %= mp_int(other);
}

inline mp_int_temp pow(const mp_int &a, int n) {
  mp_int ret;
  int prec_words = mp::prec_words;
  mp_int::mpnpwx(a, n, ret, prec_words);
  mp_int::mpnint(ret, ret, prec_words);//wastefull.
  mp_int::ovcheck(ret);
  return ret.toTempAndDestroy();
}

inline mp_int_temp pow(const mp_int &a, const mp_int& b) {
  mp_real c(pow(static_cast<mp_real>(a), static_cast<mp_real>(b)));
  int prec_words = mp::prec_words;
  // XSL if(std::abs(c[1]) > c[2] + 1.0) {
  //    mp_int::mp_int_prec_error((char*)"pow(mp_int, mp_int)");
  //  }
  mp_int d(static_cast<int>(c[0]));
  mp_real::mpnint(c, d, prec_words);
  mp_int::ovcheck(d);
  return d.toTempAndDestroy();
}

inline int integer(const mp_int& a) {
  int ret;
  ret = int(a[FST_M] * a[1]);
  //check for correctness is commented out.  Use at your own risk.
  /*  if(std::abs(a[1]) > 1.0 || ret != a) {
    cerr << (char*)"\n*** Integer conversion : A cannot be held by an integer";
    }*/
  return ret;
}

#endif 
