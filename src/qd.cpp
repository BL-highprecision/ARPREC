/*
 * src/qd.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Contains conversion routines to/from double-double and quad-double
 * numbers used in the QD library.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

#if (ARPREC_QD)
mp_real::mp_real(const dd_real &a, unsigned int sz) {
  allocate(sz);
  *this = a.x[0];
  *this += a.x[1];
}

dd_real to_dd_real(const mp_real &a) {
  double a0 = dble(a);
  double a1 = dble(a - a0);
  a0 = quick_two_sum(a0, a1, a1);
  return dd_real(a0, a1);
}

mp_real::mp_real(const qd_real &a, unsigned int sz) {
  allocate(sz);
  *this = a.x[0];
  *this += a.x[1];
  *this += a.x[2];
  *this += a.x[3];
}

qd_real to_qd_real(const mp_real &a) {
  double x0, x1, x2, x3;
  mp_real x = a;
  x0 = dble(x);
  x -= x0;
  x1 = dble(x);
  x -= x1;
  x2 = dble(x);
  x -= x2;
  x3 = dble(x);
  return qd_real(x0, x1, x2, x3);
}
#endif

