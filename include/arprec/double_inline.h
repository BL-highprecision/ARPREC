/*
 * include/arprec/double_inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2006
 *
 */
#ifndef MP_DOUBLE_INLINE_H__
#define MP_DOUBLE_INLINE_H__
/* This file contains double precision version of some function 
 * provided by the arprec library, giving compatibility when using
 * templates, etc.
 */

#include <arprec/arprec_config.h>

/* Computes the nearest integer to d. */
inline double nint(double d) {
  if (d == floor(d))
    return d;
  return floor(d + 0.5);
}

inline double aint(double a) { return a>0 ? floor(a): ceil(a); }
inline double anint(double a) { return a>0 ? ceil(a-0.5) : floor(a+0.5); }

inline void sincosh(double t, double &sinh_t, double &cosh_t) {
  sinh_t = sinh(t);
  cosh_t = cosh(t);
}

inline double sqr(double t) {
  return t * t;
}

#endif
