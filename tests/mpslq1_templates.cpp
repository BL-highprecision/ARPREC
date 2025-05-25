#ifndef MPSLQ1_TEMPLATES_CC
#define MPSLQ1_TEMPLATES_CC

#include <iomanip>
#include <stack>

#include "mpslq1.h"

template <class T>
void sort_pslq(const matrix<T> &x, matrix<int> &ip) {
  int n = x.size();
  for (int i = 0; i < n; i++)
    ip(i) = i;
  if (n <= 1)
    return;

#define SWAP(ii, jj) kk = ip(ii);  ip(ii) = ip(jj);  ip(jj) = kk;

  int kk;
  int left = 0, right = n-1;
  std::stack<int> s;
  for (;;) {
    while (right > left) {
      int m = (right + left) / 2;
      SWAP(m, left);
      int i = left - 1;
      int j = right + 1;
      const T &v = x(ip(left));
      for (;;) {
        while (x(ip(++i)) < v);
        while (x(ip(--j)) > v);
        if (i >= j) break;
        SWAP(i, j);
      }

      if (j-left > right-j) {
        s.push(left);
        s.push(j);
        left = j+1;
      } else {
        s.push(j+1); 
        s.push(right);
        right = j;
      }
    }

    if (s.empty())
      break;

    right = s.top();
    s.pop();
    left = s.top();
    s.pop();

  }
#undef SWAP

}

template <class T>
int check_history(const matrix<T> &y, const matrix<T> &hist, const T &teps) {
  int n, len;
  int result = RESULT_CONTINUE;
  T u, t;

  hist.getSize(n, len);

  for (int j = 0; j < len; j++) {
    u = 0.0;
    for (int i = 0; i < n; i++) {
      t = abs(hist(i, j) - y(i));
      if (t > u)
        u = t;
    }

    if (u < teps) {
      result = RESULT_DUPLICATE;
      if (debug_level >= 2) {
        std::cout << "Iteration " << std::setw(5) << pslq_iter << ": Duplicate found (j = " << j << ")." << std::endl;
        std::cout << "  delta = " << u << std::endl;
      }
      break;
    }
  }

  return result;
}

template <class T>
void insert_history(const matrix<T> &y, matrix<T> &hist) {
  int n, len;
  hist.getSize(n, len);
  int j = (pslq_iter % len);
  for (int i = 0; i < n; i++) {
    hist(i, j) = y(i);
  }
}

#endif

