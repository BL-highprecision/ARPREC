#ifndef PSLQ1_TEMPLATES_CC
#define PSLQ1_TEMPLATES_CC

#include <cfloat>
#include <cmath>
#include <arprec/mp_real.h>

#include "matrix.h"
#include "pslq1.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

using std::abs;

/* Compuates the smallest and largest element in magnitude 
   from the matrix v. 
   Returns the index corresponding to min if isMin is true
       otherwise the index corresponding to max
*/
template <class T>
int matrix_minmax(const matrix<T> &v, T &v_min, T &v_max, bool isMin) {
  int i, min=-1, max=-1;
  int n = v.size();

  T t;
  v_min = DBL_MAX;
  v_max = 0.0;
  for (i = 0; i < n; i++) {
    t = abs(v(i));
    if (t < v_min){
	min = i;
	v_min = t;
    }
    if (t > v_max){
	max = i;
	v_max = t;
    }
  }

  if(isMin)
      return min;
  return max;
}

/* Computes the smallest element in magnitude from matrix v. */
template <class T>
int  matrix_min(const matrix<T> &v, T &v_min) {
  int i, min=-1;
  int n = v.size();

  T t;
  v_min = DBL_MAX;
  for (i = 0; i < n; i++) {
    t = abs(v(i));
    if (t < v_min){
      v_min = t;
      min = i;
    }
  }
  return min;
}

/* Computes the largest element in magnitude from matrix v. */
template <class T>
int matrix_max(const matrix<T> &v, T &v_max) {
  int i, max;
  int n = v.size();

  T t;
  v_max = 0.0;
  for (i = 0; i < n; i++) {
    t = abs(v(i));
    if (t > v_max){
      v_max = t;
      max = i;
    }
  }
  return max;
}

/* Computes a LQ decomposition of the matrix a, and puts the lower
   triangular matrix L into a. */
template <class T>
void lq_decomp(int n, int m, matrix<T> &a) {
  //a.getSize(n, m);
  int i, j, k;
  int min_mn = std::min(m, n);
  T t, u, nrm;

  for(i = 0; i < min_mn-1; i++) {

    /* Compute the Householder vector. */
    t = a(i, i);
    nrm = sqr(t);
    for (j = i+1; j < m; j++) {
      nrm += sqr(a(i, j));
    }
    if (nrm == 0.0)
      continue;
    nrm = sqrt(nrm);
    if (t < 0.0)
      nrm = -nrm;

    t = 1.0 / nrm;
    for (j = i; j < m; j++) {
      a(i, j) *= t;
    }
    t = (a(i, i) += 1.0);

    /* Transform the rest of the rows. */
    for (j = i+1; j < n; j++) {
      u = 0.0;
      for (k = i; k < m; k++) {
        u += a(i, k) * a(j, k);
      }
      u = -u / t;

      for (k = i; k < m; k++) {
        a(j, k) += u * a(i, k);
      }
    }

    /* Set the diagonal entry.*/
    a(i, i) = -nrm;
  }

  /* Set the upper half of a to zero. */
  for (j = 0; j < m; j++) {
    for (i = 0; i < j; i++) {
      a(i, j) = 0.0;
    }
  }
}

/* Computes the bound on the relation size based on the matrix H. */
template <class T>
void bound_pslq(const matrix<T> &h, T &r) {
  int i;
  int m, n;
  h.getSize(n, m);
  int min_mn = std::min(m, n);
  T t;
  r = 0.0;
  for (i = 0; i < min_mn; i++) {
    t = abs(h(i, i));
    if (t > r)
      r = t;
  }
  r = 1.0 / r;
}

#endif
