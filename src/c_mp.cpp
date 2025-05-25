/*
 * src/c_mp.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Contains C wrapper function for ARPREC library.
 * This can be used from fortran code.
 * 
 * Last modified: July 9, 2002
 */
#include <cstdlib>
#include <string>
#include "config.h"

#include <arprec/mp_int.h>
#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include <arprec/c_mp.h>
#include <arprec/fpu.h>

#ifdef CRAY_STRINGS
#include <fortran.h>
#endif

#define DoubleToComplex(a)  (mp_complex(a, a + mp::fmpwds5))
#define MALLOC(type, n)  ((type *) std::malloc((n) * sizeof(type)))

#ifdef ALWAYS_FIX_FPU
#define FPU_FIX_START unsigned int oldcw; fpu_fix_start(&oldcw)
#define FPU_FIX_STOP  fpu_fix_end(&oldcw)
#else
#define FPU_FIX_START 
#define FPU_FIX_STOP  
#endif

using std::cout;
using std::endl;

extern "C" {

ARPREC_API int c_mpinit(int nr_digits) {
  FPU_FIX_START;
  mp::mp_init(nr_digits);
  return mp::prec_words;
  FPU_FIX_STOP;
}

ARPREC_API int c_mpinit_x(int nr_digits, const char *filename) {
  FPU_FIX_START;
  mp::mp_init(nr_digits, filename);
  return mp::prec_words;
  FPU_FIX_STOP;
}

ARPREC_API void c_mpgetpar(int pnum, int *val, int index)
{
  FPU_FIX_START;
  switch (pnum) {
  case 1: *val = mp::prec_words; break;
  case 2: *val = mp::debug_level; break;
  case 3: *val = mp::debug_words; break;
  case 5: *val = static_cast<int>(mp::round_dir); break;
  case 6: *val = mp::error_no; break;
  case 7: *val = mp::MPKER[index]; break;
  }
  FPU_FIX_STOP;
}

ARPREC_API void c_mpsetpar(int pnum, int val, int index)
{
  FPU_FIX_START;
  switch (pnum) {
  case 1: mp::prec_words = val; break;
  case 2: mp::debug_level = val; break;
  case 3: mp::debug_words = val; break;
  case 5: mp::round_dir = static_cast<enum mp::rounding_mode>(val); break;
  case 6: mp::error_no = val; break;
  case 7: mp::MPKER[index] = val; break;
  }
  FPU_FIX_STOP;
}


/* add */
ARPREC_API void c_mpadd(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  int prec_words = mp::prec_words;
  mp_real::mpadd(a2, b2, c2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_d(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) (a2 + b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_ji(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), c2(c);
  c2 = (mp_int) (a2 + b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_jd(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) (a2 + b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_zq(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5), b2(b), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpadd(ar, b2, cr, prec_words);
  ci = ai;
  ar.toTempAndDestroy();
  ai.toTempAndDestroy();
  b2.toTempAndDestroy();
  cr.toTempAndDestroy();
  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_zx(const double *a, double br, double bi, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5);
  mp_real cr(c), ci(c + mp::fmpwds5);
  cr = (mp_real) (ar + br);
  ci = (mp_real) (ai + bi);
  ar.toTempAndDestroy();
  ai.toTempAndDestroy();
  cr.toTempAndDestroy();
  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpadd_zz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5);
  mp_real br(b), bi(b + mp::fmpwds5);
  mp_real cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpadd(ar, br, cr, prec_words);
  mp_real::mpadd(ai, bi, ci, prec_words);
  ar.toTempAndDestroy(); ai.toTempAndDestroy();
  br.toTempAndDestroy(); bi.toTempAndDestroy(); 
  cr.toTempAndDestroy(); ci.toTempAndDestroy(); 
  FPU_FIX_STOP;
}

/* sub */
ARPREC_API void c_mpsub(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c); 
  int prec_words = mp::prec_words;
  mp_real::mpsub(a2, b2, c2, prec_words); 
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_d(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c); 
  c2 = (mp_real) (a2 - b); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_dq(double a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real b2(b), c2(c); 
  c2 = (mp_real) (a - b2); 
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_ji(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), c2(c);
  c2 = (mp_int) (a2 - b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_ij(int a, const double *b, double *c) {
  FPU_FIX_START;
  mp_int b2(b), c2(c);
  c2 = (mp_int) (a - b2);
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_jd(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) (a2 - b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_dj(double a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real b2(b), c2(c);
  c2 = (mp_real) (a - b2);
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_zq(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5), b2(b), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpsub(ar, b2, cr, prec_words);
  ci = ai;
  ar.toTempAndDestroy();  ai.toTempAndDestroy();
  b2.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_qz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), br(b), bi(b + mp::fmpwds5), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpsub(a2, br, cr, prec_words);
#if 0
  ci = mp_real(-bi);  //bug? -- XSL
#else
  ci = bi;
  ci[1] = -ci[1];
#endif
  a2.toTempAndDestroy();
  br.toTempAndDestroy();  bi.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_zx(const double *a, double br, double bi, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5);
  mp_real cr(c), ci(c + mp::fmpwds5);
  cr = (mp_real) (ar - br);
  ci = (mp_real) (ai - bi);
  ar.toTempAndDestroy();  ai.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_xz(double ar, double ai, const double *b, double *c) {
  FPU_FIX_START;
  mp_real br(b), bi(b + mp::fmpwds5);
  mp_real cr(c), ci(c + mp::fmpwds5);
  cr = (mp_real) (ar - br);
  ci = (mp_real) (ai - bi);
  br.toTempAndDestroy();  bi.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsub_zz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5);
  mp_real br(b), bi(b + mp::fmpwds5);
  mp_real cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpsub(ar, br, cr, prec_words);
  mp_real::mpsub(ai, bi, ci, prec_words);
  ar.toTempAndDestroy();  ai.toTempAndDestroy();
  br.toTempAndDestroy();  bi.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}


/* negation */
ARPREC_API void c_mpneg_q(const double *a, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) -a2;
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpneg_z(const double *a, double *c) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //c2 = (mp_complex) -DoubleToComplex(a);
  c2 = (mp_complex) -(a2);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* mul */
ARPREC_API void c_mpmul(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c); 
  int prec_words = mp::prec_words;
  mp_real::mpmulx(a2, b2, c2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_ji(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), c2(c);
  c2 = (mp_int) (a2 * b); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_qi(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) (a2 * (double) b); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_qd(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c); 
  int prec_words = mp::prec_words;
  mp_real::mpmuld(a2, b, 0, c2, prec_words);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_zq(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5), b2(b), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpmulx(ar, b2, cr, prec_words);
  mp_real::mpmulx(ai, b2, ci, prec_words);
  ar.toTempAndDestroy();  ai.toTempAndDestroy();
  b2.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_zz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //mp_complex::mpcmulx(DoubleToComplex(a), DoubleToComplex(b), c2);
  mp_complex::mpcmulx(a2, b2, c2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmul_zd(const double *a, double b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //mp_complex::mpcmuld(DoubleToComplex(a), b, 0, c2);
  mp_complex::mpcmuld(a2, b, 0, c2, prec_words);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* div */
ARPREC_API void c_mpdiv(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_real a2(a), b2(b), c2(c); 
  mp_real::mpdivx(a2, b2, c2, prec_words); 
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_jj(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), b2(b), c2(c);
  c2 = (mp_int) (a2 / b2); 
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_ji(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), c2(c);
  c2 = (mp_int) (a2 / b); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_ij(int a, const double *b, double *c) {
  FPU_FIX_START;
  mp_int b2(b), c2(c);
  c2 = (mp_int) ( a / b2); 
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_qi(const double *a, int b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_real a2(a), c2(c); 
  mp_real::mpdivd(a2, (double) b, 0, c2, prec_words); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_iq(int a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real b2(b), c2(c); 
  c2 = (mp_real) (a / b2); 
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_qd(const double *a, double b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_real a2(a), c2(c); 
  mp_real::mpdivd(a2, b, 0, c2, prec_words); 
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_dq(double a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real b2(b), c2(c); 
  c2 = (mp_real) (a / b2); 
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_zq(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5), b2(b), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpdivx(ar, b2, cr, prec_words);
  mp_real::mpdivx(ai, b2, ci, prec_words);
  ar.toTempAndDestroy();  ai.toTempAndDestroy();
  b2.toTempAndDestroy();
  cr.toTempAndDestroy();  ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_qz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_complex b2(b, b + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  mp_real a2(a);
  //c2 = (mp_complex) (a2 / DoubleToComplex(b));
  c2 = (mp_complex) (a2 / b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_zz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  int prec_words = mp::prec_words;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //mp_complex::mpcdivx(DoubleToComplex(a), DoubleToComplex(b), c2);
  mp_complex::mpcdivx(a2, b2, c2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_zd(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real ar(a), ai(a + mp::fmpwds5), cr(c), ci(c + mp::fmpwds5);
  int prec_words = mp::prec_words;
  mp_real::mpdivd(ar, b, 0, cr, prec_words);
  mp_real::mpdivd(ai, b, 0, ci, prec_words);
  ar.toTempAndDestroy(); ai.toTempAndDestroy();
  cr.toTempAndDestroy(); ci.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdiv_dz(double a, const double *b, double *c) {
  FPU_FIX_START;
  mp_complex b2(b, b + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //c2 = (mp_complex) (a / DoubleToComplex(b));
  c2 = (mp_complex) (a / b2);
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpdmc(double a, double *b) {
  FPU_FIX_START;
  mp_real b2(b);
  int prec_words = mp::prec_words;
  mp_real::mpdmc(a, 0, b2, prec_words);
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpmdc(const double *a, double *b, int *n) {
  FPU_FIX_START;
  mp_real a2(a);
  int prec_words = mp::prec_words;
  mp_real::mpmdc(a2, *b, *n, prec_words);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}


/*Special Functions*/
ARPREC_API void c_mperfc(const double *a, double *b) {
    mp_real a2(a), b2(b);
    b2 = erfc(a2);
    a2.toTempAndDestroy();
    b2.toTempAndDestroy();
}
ARPREC_API void c_mperf(const double *a, double *b)
{
    mp_real a2(a), b2(b);
    b2 = erf(a2);
    a2.toTempAndDestroy();
    b2.toTempAndDestroy();
}
ARPREC_API void c_mpbessel(const double *a, double *b)
{
    mp_real a2(a), b2(b);
    b2 = bessel(a2);
    a2.toTempAndDestroy();
    b2.toTempAndDestroy();
}
ARPREC_API void c_mpbesselexp(const double *a, double *b)
{
    mp_real a2(a), b2(b);
    b2 = besselexp(a2);
    a2.toTempAndDestroy();
    b2.toTempAndDestroy();
}
ARPREC_API void c_mpgamma(const double *a, double *b)
{
    mp_real a2(a), b2(b);
    b2 = gamma(a2);
    a2.toTempAndDestroy();
    b2.toTempAndDestroy();
}


/* assignment */
ARPREC_API void c_mpeq(const double *a, double *b) {
  FPU_FIX_START;
#if 0
  mp_real b2(b);
  int prec_words = mp::prec_words;
  mp_real::mpeq((mp_real) a, b2, prec_words); // cause seg. fault on Sun
  b2.toTempAndDestroy();
#else
  mp_real a2(a), b2(b);
  int prec_words = mp::prec_words;
  mp_real::mpeq(a2, b2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
#endif
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_int(int a, double *b) {
  FPU_FIX_START;
  mp_real b2(b);
  b2 = a;
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_d(double a, double *b) {
  FPU_FIX_START;
  mp_real b2(b);
  b2 = a;
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_ji(int a, double *b) {
  FPU_FIX_START;
  mp_int b2(b);
  b2 = a;
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_zq(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex b2(b, b + mp::fmpwds5);
  b2 = mp_real(a);
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_zx(double *r, double *i, double *b) {
  FPU_FIX_START;
  mp_complex b2(b, b + mp::fmpwds5);
  b2 = mp_complex(*r, *i);
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpeq_zz(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  //b2 = DoubleToComplex(a);
  b2 = a2;
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpeq_str(const char *str, double *a) {
  FPU_FIX_START;
  mp_real mp_a(a);
  mp_a = str;
  mp_a.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* power */
ARPREC_API void c_mppwr(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  c2 = (mp_real) pow((const mp_real &) a2, (const mp_real &) b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_d(const double *a, double b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  c2 = (mp_real) pow(a2, b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_qi(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  int prec_words = mp::prec_words;
  mp_real::mpnpwx(a2, b, c2, prec_words);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_jj(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), b2(b), c2(c);
  c2 = (mp_int) pow(a2, b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_ji(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_int a2(a), c2(c);
  c2 = (mp_int) pow(a2, b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_zi(const double *a, int b, double *c) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  //c2 = (mp_complex) pow(DoubleToComplex(a), b);
  c2 = (mp_complex) pow(a2, b);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mppwr_zq(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex c2(c, c + mp::fmpwds5);
  mp_real b2(b);
  // XSL ?? c2 = (mp_complex) pow(DoubleToComplex(a), b2);
  //c2 = (mp_complex) pow(DoubleToComplex(a), (mp_real) b);
  c2 = pow(a2, b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
#if 0
// bug -- XSL ??
ARPREC_API void c_mppwr_zz(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_complex c2(c, c + mp::fmpwds5);
  c2 = (mp_complex) pow(DoubleToComplex(a), DoubleToComplex(b));
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
#endif

/* equality */
ARPREC_API void c_mpcpr(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  *c = (int) (a2 == b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcpr_i(const double *a, int b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 == b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcpr_d(const double *a, double b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 == b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcpr_z(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  //c = (int) (DoubleToComplex(a) == DoubleToComplex(b));
  *c = (int) (a2 == b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* less-than-or-equal-to */
ARPREC_API void c_mplet(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  *c = (int) (a2 <= b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mplet_i(const double *a, int b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 <= b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mplet_d(const double *a, double b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 <= b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* greater-than-or-equal-to */
ARPREC_API void c_mpget(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  *c = (int) (a2 >= b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpget_i(const double *a, int b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 >= b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpget_d(const double *a, double b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 >= b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* less-than */
ARPREC_API void c_mpltt(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  *c = (int) (a2 < b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpltt_i(const double *a, int b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 < b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpltt_d(const double *a, double b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 < b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* greater-than */
ARPREC_API void c_mpgtt(const double *a, const double *b, int *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  *c = (int) (a2 > b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpgtt_i(const double *a, int b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 > b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpgtt_d(const double *a, double b, int *c) {
  FPU_FIX_START;
  mp_real a2(a);
  *c = (int) (a2 > b);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpabs(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) abs(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpabs_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_real b2(b);
  //b2 = (mp_real) abs(DoubleToComplex(a));
  b2 = (mp_real) abs(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mparg(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_real b2(b);
  //b2 = (mp_real) arg(DoubleToComplex(a));
  b2 = (mp_real) arg(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* trigonometric functions */
ARPREC_API void c_mpacos(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) acos(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpasin(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) asin(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpatan(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) atan(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpatan2(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  c2 = (mp_real) atan2(a2, b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcos(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), junk;
  mp_real::mpcssx(a2, mp_real::_pi, b2, junk);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcos_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  //b2 = (mp_complex) cos(DoubleToComplex(a));
  b2 = (mp_complex) cos(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpdble(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a);
  *b = dble(a2);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcosh(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), junk;
  mp_real::mpcshx(a2, mp_real::_pi, mp_real::_log2, b2, junk);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpexp(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  mp_real::mpexpx(a2, mp_real::_pi, mp_real::_log2, b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpexp_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  //b2 = (mp_complex) exp(DoubleToComplex(a));
  b2 = (mp_complex) exp(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpaint(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), junk(0.0, 0);
  int prec_words = mp::prec_words;
  mp_real::mpinfr(a2, b2, junk, prec_words, 0);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpnint(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  int prec_words = mp::prec_words;
  mp_real::mpnint(a2, b2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mplog(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  int prec_words = mp::prec_words;
  mp_real::mplogx(a2, mp_real::_pi, mp_real::_log2, b2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mplog_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  b2 = (mp_complex) log(a2);
  //b2 = (mp_complex) log(DoubleToComplex(a));
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mplog10(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  int prec_words = mp::prec_words;
  mp_real::mplogx(a2, mp_real::_pi, mp_real::_log2, b2, prec_words);
  mp_real::mpdivx(b2, mp_real::_log10, b2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsin(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) sin(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsin_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  b2 = (mp_complex) sin(a2);
  //b2 = (mp_complex) sin(DoubleToComplex(a));
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsinh(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) sinh(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpnrt(const double *a, int *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), c2(c);
  mp_real::mpnrtx(a2, *b, c2);
  a2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsqrt(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  int prec_words = mp::prec_words;
  mp_real::mpsqrtx(a2, b2, prec_words);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpsqrt_z(const double *a, double *b) {
  FPU_FIX_START;
  mp_complex a2(a, a + mp::fmpwds5);
  mp_complex b2(b, b + mp::fmpwds5);
  mp_complex::mpcsqrtx(a2, b2);
  //  mp_complex::mpcsqrtx(DoubleToComplex(a), b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mptan(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) tan(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mptanh(const double *a, double *b) {
  FPU_FIX_START;
  mp_real a2(a), b2(b);
  b2 = (mp_real) tanh(a2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpmod(const double *a, const double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  c2 = (mp_real) fmod(a2, b2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcsshf(const double *a, double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  mp_real::mpcshx(a2, mp_real::_pi, mp_real::_log2, b2, c2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}
ARPREC_API void c_mpcssnf(const double *a, double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  mp_real::mpcssx(a2, mp_real::_pi, b2, c2);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mprand(double *a) {
  FPU_FIX_START;
  mp_real a2(a);
  mp_real::mprand(a2);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mp_to_str(const double *a, char *str, int n_digits) {
  FPU_FIX_START;
  mp_real mp_a(a);
  std::string s = mp_a.to_string(n_digits);
  strcpy(str, s.c_str());
  mp_a.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_ovcheck(const double *a) {
  FPU_FIX_START;
  mp_int a2(a);
  mp_int::ovcheck(a2);
  a2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpinfr(const double *a, double *b, double *c) {
  FPU_FIX_START;
  mp_real a2(a), b2(b), c2(c);
  int prec_words = mp::prec_words;
  mp_real::mpinfr(a2, b2, c2, prec_words, 1);
  a2.toTempAndDestroy();
  b2.toTempAndDestroy();
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* Input */
ARPREC_API void c_mpinp(const double *q) {
  FPU_FIX_START;
  mp_real q2(q);
  std::cin >> q2;
  q2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpout(const double *q) {
  FPU_FIX_START;
  mp_real mp_q(q);
  int old_prec = cout.precision();
  cout.precision(mp::mpgetoutputprec());
  cout << mp_q;
  cout.precision(old_prec);
  mp_q.toTempAndDestroy();
  FPU_FIX_STOP;
}

/* Output */
#if 0
ARPREC_API void c_mpout(const double *q) {
  FPU_FIX_START;
  mp_real q2(q);
  cout << q2 << endl;
  q2.toTempAndDestroy();
  FPU_FIX_STOP;
}
#endif
/* c_mpwrite writes q into the string str.  On input len is the
 * maximum number of characters that should be written, including 
 * the terminating null.  On output len contains the number of 
 * characters output, not including the terminating null.   */
ARPREC_API void c_mpwrite(const double *q, char *str, int *len) {
  FPU_FIX_START;
  mp_real _q(q);
  std::string s = _q.to_string(mp::mpgetoutputprec());
  strncpy(str, s.c_str(), *len);
  str[*len-1] = 0;
  *len = static_cast<int>(s.length());
  _q.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpout_z(const double *q) {
  FPU_FIX_START;
  mp_real q2(q), q3(q + mp::fmpwds5);
  cout << "Real:" << q2 << endl;
  cout << "Imag:" << q3 << endl;
  q2.toTempAndDestroy();
  q3.toTempAndDestroy();
  FPU_FIX_STOP;
}


ARPREC_API void c_mpdotd(int *n, int *isa, double *a, int *isb, 
                         const double *db, double *c) {
  FPU_FIX_START;
  int i;
  mp_real c2(c);
#if 0
  mp_real *a2 = new mp_real[n];
#else
  mp_real *a2 = MALLOC(mp_real, *n);
#endif

#if 0
  for (i = 0; i < n; ++i) a2[i] = (mp_real) &a[i*isa];
#else
  for (i = 0; i < *n; ++i) a2[i].mpr = &a[i * *isa];
#endif

  mp_real::mpdotd(*n, 1, a2, *isb, db, c2);

#if 0
  for (i = 0; i < n; ++i) a2[i].toTempAndDestroy();
  delete [] a2;
#else
  std::free(a2);
#endif
  c2.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpsetoutputprec(int num_digits) {
  FPU_FIX_START;
  mp::mpsetoutputprec(num_digits);
  FPU_FIX_STOP;
}

ARPREC_API int c_mpgetoutputprec() {
  FPU_FIX_START;
  return mp::mpgetoutputprec();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpsetprec(int num_digits) {
  FPU_FIX_START;
  mp::mpsetprec(num_digits);
  FPU_FIX_STOP;
}

ARPREC_API int c_mpgetprec() {
  FPU_FIX_START;
  return mp::mpgetprec();
  FPU_FIX_STOP;
}

ARPREC_API void c_mpsetprecwords(int num_words) {
  FPU_FIX_START;
  mp::mpsetprecwords(num_words);
  FPU_FIX_STOP;
}

ARPREC_API int c_mpgetprecwords() {
  FPU_FIX_START;
  return mp::mpgetprecwords();
  FPU_FIX_STOP;
}

ARPREC_API void c_mppi(double *a) {
  FPU_FIX_START;
  mp_real mp_a(a);
  mp_a = mp_real::_pi;
  mp_a.toTempAndDestroy();
  FPU_FIX_STOP;
}

ARPREC_API double *c_mpnew() {
  FPU_FIX_START;
  return c_mpnew_x(mp::prec_words);
  FPU_FIX_STOP;
}

ARPREC_API double *c_mpnew_x(int nr_words) {
  FPU_FIX_START;
  int n_words = nr_words + 5;
  double *a = MALLOC(double, n_words);
  a[0] = (double) n_words;
  a[1] = (double) nr_words;
  return a;
  FPU_FIX_STOP;
}

ARPREC_API void c_mpfree(double *a) {
  FPU_FIX_START;
  std::free(a);
  FPU_FIX_STOP;
}

}  // end extern
