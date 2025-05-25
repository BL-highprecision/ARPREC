#ifndef PSLQ2_TEMPLATES_CC
#define PSLQ2_TEMPLATES_CC

#include <iostream>
#include "pslq2.h"

using std::abs;

/* Multiples a m-by-n matrix B to the left by a square matrix A. 
   The result is put into B. */
template <class T>
void matmul_left(const matrix<T> &a, matrix<mp_real> &b) {
  int m, n;
  int i, j;

  b.getSize(m, n);
  matrix<mp_real> c(m, n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      pslq_dot(m, 1, &b(0, j), m, &a(i), c(i, j));
    }
  }

  b = c;
}

/* Multiplies m-by-n matrix a by square matrix b, putting the result into a. */
template <class T>
void matmul_right(matrix<mp_real> &a, const matrix<T> &b) {
  int m, n;
  int i, j;

  a.getSize(m, n);
  matrix<mp_real> c(m, n);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      pslq_dot(n, m, &a(i), 1, &b(0, j), c(i, j));
    }
  }

  a = c;
}

/* Multiplies the transpose of the m-by-n matrix a by square matrix b, 
   putting the result into transpose of a. */
template <class T>
void matmul_right_trans(matrix<mp_real> &a, const matrix<T> &b) {
  int m, n;
  int i, j;

  a.getSize(m, n);
  matrix<mp_real> c(m, n);

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      pslq_dot(m, 1, &a(0, i), 1, &b(0, j), c(j, i));
    }
  }

  a = c;
}

/* Update the higher precision matrices B, H, and vector y using 
   a lower precision matrices dA, dB, dH, and dy.

   Returns RESULT_RELATION_FOUND is the minimum element of the
   updated matrix drops below the eps.  Returns RESULT_PRECISION_EXHAUSTED
   if the minimum is below teps but above eps.  */
template <class T>
int update_pslq(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &b, 
                matrix<mp_real> &h, matrix<mp_real> &y, const mp_real &eps, 
                const mp_real &teps) {
  int result = RESULT_CONTINUE;
  mp_real u, v;

  /* Compute y = y * db */
  matmul_right_trans(y, db);

  matrix_minmax(y, u, v);
  if (u < teps) {
    result = (u < eps) ? RESULT_RELATION_FOUND : RESULT_PRECISION_EXHAUSTED;
  }

  /* Compute b = b * db */
  matmul_right(b, db);

  /* Compute h = da * h */
  matmul_left(da, h);

  return result;
}

template <class T>
int update_pslq2(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &a, 
                 matrix<mp_real> &b, matrix<mp_real> &h, matrix<mp_real> &y, 
                 const mp_real &eps, const mp_real &teps) {
  int result = update_pslq(da, db, b, h, y, eps, teps);
  matmul_left(da, a);
  return result;
}

template <class T>
int check_pslq(const matrix<mp_real> &y, const T &eps, int nr_words) {
  PREC_START;

  int result = RESULT_CONTINUE;
  mp_real u, v, t;

  matrix_minmax(y, u, v);

  t = u / v;
  if (t < eps) {
    result = RESULT_PRECISION_EXHAUSTED;
  }

  PREC_END;
  return result;
}

/* Performs one iteration of PSLQ in double precision arithmetic. */
template <class T>
int iterate_pslq2(double gamma, matrix<T> &a, matrix<T> &b, 
                  matrix<T> &h, matrix<T> &y, 
                  const T &teps, const T &mx1, const T &mx2, int nr_words) { 
  PREC_START;
  int i, j, k;
  int n = y.size();
  T t1, t2, t3, t4;
  double d;
  int im, im1;

  /* Find the diagonal element h(j, j) such that 
     |gamma^j h(j, j)| is maximized.                   */
  t1 = 0.0;
  im = -1;
  d = gamma;
  for (i = 0; i < n-1; i++, d *= gamma) {
    t2 = d * abs(h(i, i));
    if (t2 > t1) {
      im = i;
      t1 = t2;
    }
  }

  if (im == -1) {
    std::cerr << "ERROR: Invalid index." << std::endl;
    exit(-1);
  }

  /* Exchange the im and im+1 entries of y, rows of h, 
     and columns of a and b. */
  im1 = im + 1;

  t1 = y(im);
  y(im) = y(im1);
  y(im1) = t1;

  for (i = 0; i < n-1; i++) {
    t1 = h(im, i);
    h(im, i) = h(im1, i);
    h(im1, i) = t1;
  }

  for (i = 0; i < n; i++) {
    t1 = a(im, i);
    a(im, i) = a(im1, i);
    a(im1, i) = t1;
    t1 = b(i, im);
    b(i, im) = b(i, im1);
    b(i, im1) = t1;
  }

  /* Update H with permutation produced above. */
  if (im <= n-3) {
    t1 = h(im, im);
    t2 = h(im, im1);
    t3 = 1.0 / sqrt(sqr(t1) + sqr(t2));
    t1 *= t3;
    t2 *= t3;

    for (i = im; i < n; i++) {
      t3 = h(i, im);
      t4 = h(i, im1);
      h(i, im) = t1 * t3 + t2 * t4;
      h(i, im1) = t1 * t4 - t2 * t3;
    }
  }

  /* Reduce H, updating y, A, B, and H.  We also keep track of 
     largest entry in A and B. */
  t2 = 0.0;
  for (i = im1; i < n; i++) {
    int j1 = (i == im1) ? i-1 : im1;
    
    for (j = j1; j >= 0; j--) {
      t1 = anint(h(i, j) / h(j, j));
      if (t1 == 0.0)
        continue;

      y(j) += t1 * y(i);

      for (k = 0; k < n; k++) {
        a(i, k) -= t1 * a(j, k);
        if (a(i, k) > t2)
          t2 = a(i, k);

        b(k, j) += t1 * b(k, i);
        if (b(k, j) > t2)
          t2 = b(k, j);
      }

      for (k = 0; k <= j; k++) {
        h(i, k) -= t1 * h(j, k);
      }
    }
  }

  /* Find the min of |y|. */
  matrix_min(y, t1);

  int result = RESULT_CONTINUE;
  if (t1 < teps) {
    result = RESULT_RELATION_FOUND;
  }
  
  if (t2 > mx1) {
    if (t2 > mx2) {
      result = RESULT_VERY_LARGE_VALUE;
    } else {
      result = RESULT_LARGE_VALUE;
    }
  }

  PREC_END;
  return result;
}

inline void pslq_round(mp_real &s, const mp_real &a) {
  s = a;
}

inline void pslq_round(double &s, const mp_real &a) {
  s = dble(a);
}

template <class T>
int init_pslq2(matrix<T> &dy, matrix<T> &dh, matrix<T> &da, matrix<T> &db, 
               const matrix<mp_real> &y, const matrix<mp_real> &h, 
               const T &teps, int nr_words) {
  PREC_START;

  int i, j;
  int n = y.size();
  int result = RESULT_CONTINUE;
  mp_real u, v, t;
  
  /* Find the min and max magnitude in the y vector.  If the dynamic
     range of the y vector is too large, abort the initialization.    */
  matrix_minmax(y, u, v);
  t = u / v;
  if (t < teps) {
    result = RESULT_RANGE_TOO_LARGE;
  } else {

    /* Set dy to be the scaled y vector. */
    t = 1.0 / v;
    for (i = 0; i < n; i++) {
      pslq_round(dy(i), mp_real(y(i) * t));
    }

    /* Find the max magnitude along diagonal of H. */
    v = 0.0;
    for (i = 0; i < n-1; i++) {
      t = abs(h(i, i));
      if (t > v)
        v = t;
    }

    /* Set dh to be the scaled H matrix. */
    t = 1.0 / v;
    for (j = 0; j < n-1; j++) {
      for (i = 0; i < n; i++) {
        pslq_round(dh(i, j), mp_real(h(i, j) * t));
      }
    }

    /* Set da and db to the identity. */
    da.identity();
    db.identity();
  }

  PREC_END;
  return result;
}

/* Copies the matrices a, b, h and y into sa, sb, sh, and sy, 
   respectively.  The matrices a, b, and h are n-by-n, while
   y is a vector of length n.  */
template <class T>
void save_pslq(const matrix<T> &a, const matrix<T> &b, 
               const matrix<T> &h, const matrix<T> &y, 
               matrix<T> &sa, matrix<T> &sb, matrix<T> &sh, 
               matrix<T> &sy) {
  sa = a;
  sb = b;
  sh = h;
  sy = y;
}

#endif

