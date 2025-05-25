/*
 * include/arprec/mp_complex_temp.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MP_COMPLEX_TEMP_H
#define __MP_COMPLEX_TEMP_H

class mp_complex_temp: public mp
{
 protected:
 public:
  mp_real_temp real, imag;
  mp_complex_temp(int s = n_words): real(s), imag(s) {}
  mp_complex_temp(const mp_real_temp& a, const mp_real_temp& b)
    : real(a), imag(b) {}
  mp_complex_temp(double *r, double *i): real(r), imag(i) {}
  mp_complex_temp(const double *r, const  double *i): real(r), imag(i) {}
};


#endif /* __MP_COMPLEX_TEMP_H */
