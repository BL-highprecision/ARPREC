#ifndef PSLQ2_H
#define PSLQ2_H

#include "matrix.h"
#include "pslq1.h"

/* Level 2 routines.  Found in pslq_utils2.cpp. */
template <class T>
int check_pslq(const matrix<mp_real> &y, const T &eps, int nr_words = 0);

template <class T>
int update_pslq(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &b, 
                matrix<mp_real> &h, matrix<mp_real> &y, const mp_real &eps, 
                const mp_real &teps);
template <class T>
int update_pslq2(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &a, 
                 matrix<mp_real> &b, matrix<mp_real> &h, matrix<mp_real> &y, 
                 const mp_real &eps, const mp_real &teps);

template <class T>
int iterate_pslq2(double gamma, matrix<T> &a, matrix<T> &b, matrix<T> &h, 
                  matrix<T> &y, const T &teps, const T &mx1, const T &mx2, 
                  int nw_words = 0);

inline int iterate_pslq2(double gamma, matrix<double> &a, 
                  matrix<double> &b, matrix<double> &h, matrix<double> &y) {
  return iterate_pslq2(gamma, a, b, h, y, 
      1.0e-14, 1.0e13, 4.503599627370496e15);
}

template <class T>
int init_pslq2(matrix<T> &dy, matrix<T> &dh, matrix<T> &da, 
               matrix<T> &db, const matrix<mp_real> &y, 
               const matrix<mp_real> &h, const T &teps, int nr_words = 0);

template <class T>
void save_pslq(const matrix<T> &a, const matrix<T> &b, 
               const matrix<T> &h, const matrix<T> &y, 
               matrix<T> &sa, matrix<T> &sb, matrix<T> &sh, 
               matrix<T> &sy);

int pslq2(const matrix<mp_real> &x, matrix<mp_real> &rel, 
          const mp_real &eps, double gamma = DEFAULT_GAMMA);

void pslq_dot(int n, int inca, const mp_real *a, 
              int incb, const double *b, mp_real &c);
void pslq_dot(int n, int inca, const mp_real *a, 
              int incb, const mp_real *b, mp_real &c);

template <class T> 
void matmul_right(matrix<mp_real> &a, matrix<T> &b);
template <class T> 
void matmul_left(matrix<T> &a, matrix<mp_real> &b);
template <class T> 
void matmul_left_trans(matrix<T> &a, matrix<mp_real> &b);

#include "pslq2_templates.cpp"

#endif
