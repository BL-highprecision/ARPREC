#ifndef MPSLQ2_H_
#define MPSLQ2_H_

#include "pslq2.h"
#include "mpslq1.h"

template <class T>
int iterate_mpslq2(double gamma, matrix<T> &a, matrix<T> &b, 
		               matrix<T> &h, matrix<T> &y, const T &teps, 
                   bool use_only_one_pair, 
									 const T &mx1, const T &mx2, int nr_words = 0);

inline int iterate_mpslq2(double gamma, matrix<double> &a, 
                          matrix<double> &b, matrix<double> &h, 
													matrix<double> &y, bool use_only_one_pair) {
  return iterate_mpslq2(gamma, a, b, h, y, 
      1.0e-14, use_only_one_pair, 1.0e13, 4.503599627370496e15);
}

template <class T>
int update_mpslq(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &b, 
                 matrix<mp_real> &h, matrix<mp_real> &y, const mp_real &eps, 
                 const mp_real &teps);

int mpslq2(const matrix<mp_real> &x, matrix<mp_real> &rel, 
           const mp_real &eps, double gamma = DEFAULT_GAMMA);

#include "mpslq2_templates.cpp"
#endif

