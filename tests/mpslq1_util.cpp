#include "mpslq1.h"

using std::cerr;
using std::cout;
using std::endl;

int iterate_mpslq(double gamma, matrix<mp_real> &y, matrix<mp_real> &h, 
                  matrix<mp_real> &b, const mp_real &eps, const mp_real &teps, 
                  bool use_only_one_pair) {
  int i, j;
  int n = y.size();
  double beta = 0.4;

  int mnp = use_only_one_pair ? 1 : (int) (beta * n);
  if (mnp < 1)
    mnp = 1;

  /* Compute and sort abs(gamma ^ i * h(i, i)) in increasing order. */
  matrix<mp_real> q(n-1);
  matrix<int> ip(n-1);
  matrix<bool> ib(n);
  matrix<int> is(mnp);

  double d = gamma;
  for (i = 0; i < n-1; i++, d *= gamma) {
    q(i) = abs(d * h(i, i));
  }
  sort_pslq(q, ip);

  /* Pick up to mnp pairs (im, im+1), starting with largest. */
  for (i = 0; i < n; i++)
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
  mp_real u;
  for (i = 0; i < np; i++) {
    int im = is(i);
    int im1 = im + 1;

    u = y(im);
    y(im) = y(im1);
    y(im1) = u;

    for (j = 0; j < n; j++) {
      u = b(im, j);
      b(im, j) = b(im1, j);
      b(im1, j) = u;
    }

    for (j = 0; j < n-1; j++) {
      u = h(im, j);
      h(im, j) = h(im1, j);
      h(im1, j) = u;
    }
  }

  /* Remove corners */
  mp_real t1, t2, t3, t4;
  for (i = 0; i < np; i++) {
    int im = is(i);
    int im1 = im + 1;
    if (im == n-2)
      continue;

    t1 = h(im, im);
    t2 = h(im, im1);
    t3 = 1.0 / sqrt(sqr(t1) + sqr(t2));
    t1 *= t3;
    t2 *= t3;

    for (j = im; j < n; j++) {
      t3 = h(j, im);
      t4 = h(j, im1);
      h(j, im) = t1 * t3 + t2 * t4;
      h(j, im1) = t1 * t4 - t2 * t3;
    }
  }

  /* Reduce H. */
  int r;
  matrix<mp_real> t(n, n);
  for (i = 1; i < n; i++) {
    for (j = 0; j < n-i; j++) {
      int ij = i + j;
      for (r = j+1; r < ij; r++) {
        h(ij, j) -= t(ij, r) * h(r, j);
      }

      t(ij, j) = 
        anint(h(ij, j) / h(j, j));
      h(ij, j) -= t(ij, j) * h(j, j);
    }
  }

  /* Update y, using the T array.  The smallest magnitude element of
     y is found at the same time. */
  t1 = abs(y(n-1));
  int ymin = -1;
  for (j = 0; j < n-1; j++) {
    for (i = j+1; i < n; i++) {
      y(j) += t(i, j) * y(i);
    }
    t2 = abs(y(j));
    if (t2 < t1){
      t1 = t2;
      ymin = j;
    }
  }

  if(ymin < 0){
      cerr << "mpslq1_util.cpp: Error ymin=-1" << endl;
      exit(-1);
  }

  /* Update b, using the T array. */
  for (r = 0; r < n; r++) {
    for (j = 0; j < n-1; j++) {
      for (i = j + 1; i < n; i++) {
        b(j, r) += t(i, j) * b(i, r);
      }
    }
  }

  t2 = 0.0;
  for(int i =0; i< n; i++){
      t3 = abs(b(ymin,i));
      t2 = std::max(t2,t3);
  }

  int result = RESULT_CONTINUE;
  if (t1 < t2*teps) {
    result = (t1 < t2*eps) ? RESULT_RELATION_FOUND : RESULT_PRECISION_EXHAUSTED;
  }

  return result;
}

