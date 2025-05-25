#include <iostream>
#include <iomanip>
#include <cmath>

#include <arprec/mp_real.h>

#include "pslq1.h"

using std::cerr;
using std::cout;
using std::endl;

int pslq1(const matrix<mp_real> &x, matrix<mp_real> &rel, 
          const mp_real &eps, double gamma) {
  int print_interval = 100;
  int check_interval = 500;
  int n = x.size();
  matrix<mp_real> b(n, n);
  matrix<mp_real> h(n, n-1);
  matrix<mp_real> y(n);
  mp_real t;
  mp_real teps = eps * 1.0e20;
  mp_real max_bound = 0.0;
  std::ios_base::fmtflags fmt = cout.flags();
  cout << std::scientific << std::setprecision(20);

  int result;

  init_pslq(x, y, h, b);
  if (debug_level >= 3) {
    y.print("Initial y:");
    b.print("Initial B:");
    h.print("Initial H:");
  }
  result = reduce_pslq(h, y, b, eps);

  pslq_iter = 0;
  while (result == RESULT_CONTINUE) {
    pslq_iter++;
    if (debug_level >= 2) {
      if (pslq_iter % print_interval == 0)
        cout << "Iteration " << std::setw(5) << pslq_iter << endl;
    }

    result = iterate_pslq(gamma, y, h, b, eps, teps);

    if (debug_level >= 3) {
      y.print("Updated y: ");
      b.print("Updated B: ");
      h.print("Updated H: ");
    }

    if (pslq_iter % check_interval == 0) {
      /* Find min and max magnitude in y vector. */
      mp_real u, v;
      matrix_minmax(y, u, v);

      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(5) << pslq_iter << endl;
        cout << "  min(y) = " << u << endl;
        cout << "  max(y) = " << v << endl;
      }

      /* Compute norm bound. */
      bound_pslq(h, u);
      if (u > max_bound)
        max_bound = u;

      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(5) << pslq_iter << endl;
        cout << "  norm bound = " << u << endl;
        cout << "  max  bound = " << max_bound << endl;
      }
    }

  } 

  int jm = -1;
  if (result == RESULT_RELATION_FOUND) {
    mp_real u, v;

    /*Output final norm bound.*/
    /*Relation found.  Select the relation with smallest y and compute norm.*/
    jm = matrix_minmax(y, u, v);
    bound_pslq(h, t);
    cout << "Relation detected at iteration " << pslq_iter << endl;
    cout << "  min(y) = " << u << endl;
    cout << "  max(y) = " << v << endl;
    cout << "  bound  = " << t << endl;

    if (jm < 0) {
      cerr <<  "ERROR: Invalid index." << endl;
      exit(-1);
    }

    for (int i = 0; i < n; i++) {
      rel(i) = b(i, jm);
    }
  }

  cout.flags(fmt);
  return result;
}

