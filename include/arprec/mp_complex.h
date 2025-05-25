/*
 * include/arprec/mp_complex.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MPCOMPLEX_H
#define __MPCOMPLEX_H

#include <arprec/mp_real.h>
#include <arprec/mp_complex_temp.h>

struct ARPREC_API mp_complex : mp {
  mp_real real, imag;

  mp_complex():real(), imag() {}
  mp_complex(const mp_real& a, const mp_real& b):real(a), imag(b) {}
  mp_complex(const mp_real& a) : real(a) { mp_real::zero(imag); }
  mp_complex(const mp_complex& x) : real(x.real), imag(x.imag) {}
  mp_complex(double dpr, double dpi = 0.0, int s = n_words)
            : real(dpr, s), imag(dpi, s) {}
  mp_complex(double *r, double *i) : real(r), imag(i) {}
  mp_complex(const double *r, const double *i) : real(r), imag(i) {}
  mp_complex(const mp_complex_temp& x) : real(x.real), imag(x.imag) {}

  mp_complex_temp toTempAndDestroy() { 
    return mp_complex_temp(real.toTempAndDestroy(), imag.toTempAndDestroy());
  }
  
  static void zero(mp_complex& a) { 
    mp_real::zero(a.real); 
    mp_real::zero(a.imag); 
  }

  static void mpcmuld(const mp_complex& a, double db, int n, 
		      mp_complex& b, int prec_words);
  static void mpcadd(const mp_complex& a, const mp_complex& b, mp_complex& c,
		     int prec_words);
  static void mpcsub(const mp_complex& a, const mp_complex& b, mp_complex& c, 
		     int prec_words);
  static void mpcmul(const mp_complex& a, const mp_complex& b, mp_complex& c,
		     int prec_words);
  static void mpcmulx(const mp_complex& a, const mp_complex& b, mp_complex& c,
		      int prec_words);
  static void mpcdiv(const mp_complex& a, const mp_complex& b,
		     mp_complex& c, int prec_words);
  static void mpcdivx(const mp_complex& a, const mp_complex& b,
		      mp_complex& c, int prec_words);
  static void mpceq(const mp_complex& a, mp_complex& b, int prec_words);
  static void mpcsqx(const mp_complex& a, mp_complex& c, int prec_words);
  static void mpcagx(mp_complex& a, mp_complex& b);
  static void mpcsqrtx(const mp_complex& a, mp_complex& b);
  static void mpcpwx(const mp_complex& a, int n, mp_complex& b);
  static void mpcpwr(const mp_complex& a, int n, mp_complex& b);
  static void mpcsqrt(const mp_complex &a, mp_complex& b);


  mp_complex &operator+=(const mp_complex &x);
  mp_complex &operator-=(const mp_complex &x);
  mp_complex &operator*=(const mp_complex &x);
  mp_complex &operator/=(const mp_complex &x);

  mp_complex &operator+=(const mp_real &x);
  mp_complex &operator-=(const mp_real &x);
  mp_complex &operator*=(const mp_real &x);
  mp_complex &operator/=(const mp_real &x);

  mp_complex &operator/=(double x);
  mp_complex &operator/=(int x);

  mp_complex &operator=(const mp_complex &x);
  mp_complex &operator=(const mp_real &x);
  mp_complex &operator=(mp_complex_temp &x);

  mp_complex_temp operator-() {
    return mp_complex_temp(-(this->real), -(this->imag));
  }

  static mp_complex *alloc_array(int n);
  static void free_array(mp_complex *m);

};

ARPREC_API mp_complex_temp operator+(const mp_complex &a, 
          const mp_complex &b);
ARPREC_API mp_complex_temp operator+(const mp_complex &a, const mp_real &breal);
ARPREC_API mp_complex_temp operator+(const mp_real &breal, const mp_complex &a);

ARPREC_API mp_complex_temp operator-(const mp_complex &a, 
         const mp_complex &b);
ARPREC_API mp_complex_temp operator-(const mp_complex &a, 
         const mp_real &breal);
ARPREC_API mp_complex_temp operator-(const mp_real &breal, 
         const mp_complex &a);

ARPREC_API mp_complex_temp operator*(const mp_complex &a, 
          const mp_complex &b);
ARPREC_API mp_complex_temp operator/(const mp_complex &a, 
         const mp_complex &b);
ARPREC_API mp_complex_temp operator/(const mp_complex &a, 
         const mp_real &b);
ARPREC_API mp_complex_temp operator/(const mp_real &a, 
         const mp_complex &b);
ARPREC_API mp_complex_temp operator/(const mp_complex &a, 
         double b);
ARPREC_API mp_complex_temp operator/(double b, const mp_complex &a);
ARPREC_API mp_complex_temp operator/(int b, const mp_complex &a);

ARPREC_API mp_complex_temp operator*(const mp_complex& a, const mp_real& b);
ARPREC_API mp_complex_temp operator*(const mp_real& b, const mp_complex& a);
ARPREC_API mp_complex_temp operator*(const mp_complex& a, double b);
ARPREC_API mp_complex_temp operator*(double b, const mp_complex& a);
ARPREC_API mp_complex_temp exp(const mp_complex& a);
ARPREC_API mp_complex_temp log(const mp_complex& a);
ARPREC_API mp_complex_temp sin(const mp_complex &a);
ARPREC_API mp_complex_temp cos(const mp_complex &a);
ARPREC_API mp_complex_temp sqr(const mp_complex &a);
ARPREC_API mp_complex_temp sqrt(const mp_complex &a);
ARPREC_API mp_real_temp abs(const mp_complex &a);
ARPREC_API mp_real_temp arg(const mp_complex &a);
ARPREC_API mp_complex_temp pow(const mp_complex& a, int n);
ARPREC_API mp_complex_temp pow(const mp_complex& a, const mp_real& b);
ARPREC_API mp_complex_temp pow(const mp_complex& a, const mp_complex& b);

ARPREC_API bool operator==(const mp_complex &a, const mp_complex &b);
ARPREC_API bool operator!=(const mp_complex &a, const mp_complex &b);


#if (ARPREC_INLINE)
#include <arprec/mp_complex_inline.h>
#else
mp_complex_temp operator+(const mp_complex &a, const mp_complex &b);
mp_complex_temp operator+(const mp_complex &a, const mp_real &breal);
mp_complex_temp operator+(const mp_real &breal, const mp_complex &a);
mp_complex_temp operator-(const mp_complex &a, const mp_complex &b);
mp_complex_temp operator-(const mp_complex &a, const mp_real &breal);
mp_complex_temp operator-(const mp_real &breal, const mp_complex &a);
mp_complex_temp operator*(const mp_complex &a, const mp_complex &b);
mp_complex_temp operator/(const mp_complex &a, const mp_complex &b);
mp_complex_temp operator/(const mp_complex &a, const mp_real &b);
mp_complex_temp operator/(const mp_complex &a, double b);
mp_complex_temp operator/(const mp_real &a, const mp_complex &b);
mp_complex_temp operator/(double b, const mp_complex& a);
mp_complex_temp operator/(int b, const mp_complex& a);
mp_complex_temp operator*(const mp_complex& a, double b);
mp_complex_temp operator*(double b, const mp_complex& a);
mp_complex_temp operator*(const mp_complex& a, const mp_real& b);
mp_complex_temp operator*(const mp_real& b, const mp_complex& a);
mp_complex_temp exp(const mp_complex& a);
mp_complex_temp log(const mp_complex& a);
mp_complex_temp sin(const mp_complex &a);
mp_complex_temp cos(const mp_complex &a);
mp_complex_temp cos(const mp_complex &a);
bool operator==(const mp_complex &a, const mp_complex &b);
bool operator!=(const mp_complex &a, const mp_complex &b);
mp_complex_temp sqr(const mp_complex &a);
mp_complex_temp sqrt(const mp_complex &a);
mp_real_temp abs(const mp_complex &a);
mp_real_temp arg(const mp_complex &a);
mp_complex_temp pow(const mp_complex& a, int n);
mp_complex_temp pow(const mp_complex& a, const mp_real& b);
mp_complex_temp pow(const mp_complex& a, const mp_complex& b);

#endif

#endif /* __MPCOMPLEX_H */
