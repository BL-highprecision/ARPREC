/*
 * include/arprec/mp_int_temp.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MPINT_TEMP_H
#define __MPINT_TEMP_H

#include <arprec/mp_real_temp.h>
#include <arprec/mp_real.h>

class mp_int_temp: public mp
{
 protected:
 public:
  double *mpr;
  mp_int_temp(int s = n_words) {
    mpr = new double[s];
    mpr[0] = s;
  }
  mp_int_temp(double *dp) { mpr = dp; }
  mp_int_temp(const double *dp) { mpr = const_cast<double *>(dp); }
  mp_int_temp(const mp_int_temp& other) { mpr = other.mpr; } 
  mp_int_temp(const mp_real_temp& other) { mpr = other.mpr; } 

  const mp_real_temp to_mp_real_temp() const { return mp_real_temp(mpr); }
  // intentionally no deconstructor -- might use reference counting later.

};

#endif
