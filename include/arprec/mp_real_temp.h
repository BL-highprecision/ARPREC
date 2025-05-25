/*
 * include/arprec/mp_real_temp.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MPREAL_TEMP_H
#define __MPREAL_TEMP_H

#include <cstddef>

class mp_real_temp
{
 protected:
 public:
  double *mpr;
  mp_real_temp(int s = mp::n_words) {
    mpr = new double[s];
    mpr[0] = s;
  }
  mp_real_temp(double *dp) { mpr = dp; }
  mp_real_temp(const double *dp) { mpr = const_cast<double *>(dp); }
  mp_real_temp(const mp_real_temp& other) { mpr = other.mpr; } 
  // intentionally no deconstructor -- might use reference counting later.

  //inline operator double() const;
};

#endif
