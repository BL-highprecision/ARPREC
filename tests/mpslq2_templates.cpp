#ifndef MPSLQ2_TEMPLATES_CC
#define MPSLQ2_TEMPLATES_CC

#include <arprec/mp_real.h>
#include "matrix.h"
#include "mpslq2.h"

template <class T>
int iterate_mpslq2(double gamma, matrix<T> &a, matrix<T> &b, matrix<T> &h, 
                   matrix<T> &y, const T &teps, bool use_only_one_pair, 
									 const T &mx1, const T &mx2, int nr_words) {
	PREC_START;
	int n = y.size();
  double beta = 0.4;

  int mnp = use_only_one_pair ? 1 : (int) (beta * n);
  if (mnp < 1)
    mnp = 1;

  /* Compute and sort abs(gamma ^ i * h(i, i)) in increasing order. */
  matrix<T> q(n-1);
	matrix<int> ip(n-1);
	matrix<bool> ib(n);
	matrix<int> is(mnp);

  double d = gamma;
  for (int i = 0; i < n-1; i++, d *= gamma) {
    q(i) = abs(d * h(i, i));
  }
  sort_pslq(q, ip);

  /* Pick up to mnp pairs (im, im+1), starting with largest. */
  for (int i = 0; i < n; i++)
    ib(i) = false;
  int np = 0;
  int k = n-2;
  while (np < mnp) {
    int kk = ip(k);
    if (ib(kk) || ib(kk+1)) {
      if (--k < 0)
        break;
      continue;
    }
    is(np++) = kk;
    ib(kk) = ib(kk+1) = true;
  }

  /* Exchange rows. */
  T u;
  for (int i = 0; i < np; i++) {
    int im = is(i);
    int im1 = im + 1;

    u = y(im);
    y(im) = y(im1);
    y(im1) = u;

    for (int j = 0; j < n; j++) {
      u = b(im, j);
      b(im, j) = b(im1, j);
      b(im1, j) = u;

      u = a(im, j);
      a(im, j) = a(im1, j);
      a(im1, j) = u;
    }

    for (int j = 0; j < n-1; j++) {
      u = h(im, j);
      h(im, j) = h(im1, j);
      h(im1, j) = u;
    }
  }

  /* Remove corners */
  T t1, t2, t3, t4;
  for (int i = 0; i < np; i++) {
    int im = is(i);
    int im1 = im + 1;
    if (im == n-2)
      continue;

    t1 = h(im, im);
    t2 = h(im, im1);
    t3 = 1.0 / sqrt(sqr(t1) + sqr(t2));
    t1 *= t3;
    t2 *= t3;

    for (int j = im; j < n; j++) {
      t3 = h(j, im);
      t4 = h(j, im1);
      h(j, im) = t1 * t3 + t2 * t4;
      h(j, im1) = t1 * t4 - t2 * t3;
    }
  }

  /* Reduce H. */
  matrix<T> t(n, n);
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < n-i; j++) {
      int ij = i + j;
      for (int k = j+1; k < ij; k++) {
        h(ij, j) -= t(ij, k) * h(k, j);
      }

      t(ij, j) = 
        anint(h(ij, j) / h(j, j));
      h(ij, j) -= t(ij, j) * h(j, j);
    }
  }

  /* Update y, using the T array.  The smallest magnitude element of
     y is found at the same time. */
  t1 = abs(y(n-1));
  for (int j = 0; j < n-1; j++) {
    for (int i = j+1; i < n; i++) {
      y(j) += t(i, j) * y(i);
    }
    t2 = abs(y(j));
    if (t2 < t1)
      t1 = t2;
  }

  /* Update b, using the T array. */
	t2 = 0.0;
  for (int k = 0; k < n; k++) {
    for (int j = 0; j < n-1; j++) {
      for (int i = j + 1; i < n; i++) {
        b(j, k) += t(i, j) * b(i, k);
	a(i, k) -= t(i, j) * a(j, k);
	t3 = abs(b(j, k));
	if (t3 > t2)
	    t2 = t3;
	t3 = abs(a(i, k));
	if (t3 > t2)
	    t2 = t3;
      }
    }
  }

  int result = RESULT_CONTINUE;
  if (t1 < teps) {
    result = RESULT_RELATION_FOUND;
  }

  if (t2 > mx1) {
      result = (t2 > mx2) ? RESULT_VERY_LARGE_VALUE : RESULT_LARGE_VALUE;
  }

  PREC_END;
  return result;
}

template <class T>
int update_mpslq(const matrix<T> &da, const matrix<T> &db, matrix<mp_real> &b, 
                 matrix<mp_real> &h, matrix<mp_real> &y, const mp_real &eps, 
                 const mp_real &teps) {
  int result = RESULT_CONTINUE;
  mp_real u, v;

  /* Compute y = y * db */
  matmul_left(db, y);

    /* Compute b = b * db */
  matmul_left(db, b);

  /* Compute h = da * h */
  matmul_left(da, h);

  int jm = matrix_minmax(y, u, v), size=y.size();
  mp_real bMax = 0.0, t;
  for(int i=0; i<size; i++){
      t = abs(b(jm,i));
      bMax = std::max(bMax, t);
  }

  if (u < bMax*teps) {
    result = (u < bMax*eps) ? RESULT_RELATION_FOUND : RESULT_PRECISION_EXHAUSTED;
  }

  return result;
}


#endif
