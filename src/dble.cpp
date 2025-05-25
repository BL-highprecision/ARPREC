/*
 * src/dble.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002-2008
 *
 */
#include <arprec/mp_real.h>
#include <limits>

static const double _d_inf = std::numeric_limits<double>::infinity();

double dble(const mp_real &a) {
  if (a[1] == 0.0) { return 0.0; }
  if (a[2] >= 22.0) { return (a[1] > 0.0) ? _d_inf : -_d_inf; }
  if (a[2] <= -24.0) { return (a[1] > 0.0) ? 0.0 : -0.0; }

  int na = static_cast<int>(std::abs(a[1]));
  double da = a[FST_M];
  if (na == 2) {
    da += a[FST_M + 1] * mp::mprdx;
  } else if (na >= 3) {
    // we need to add from smallest to largest
    da += (a[FST_M + 1] * mp::mprdx + a[FST_M + 2] * mp::mprx2);
  }
  
  if (a[2] == -23.0) {
    // We are close to double precision underflow; make sure 
    // ldexp does not underflow to zero.
    da *= mp::mprdx;
    da *= std::ldexp(1.0, -mp::mpnbt * 22);
  } else {
    da *= std::ldexp(1.0, mp::mpnbt * static_cast<int>(a[2]));
  }
  return (a[1] > 0.0) ? da : -da;
}

