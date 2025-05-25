#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <arprec/mp_real.h>
#include <arprec/mp_int.h>

#include "mpslq1.h"

using std::cerr;
using std::cout;
using std::endl;

int mpslq1(const matrix<mp_real> &x, matrix<mp_real> &rel, 
           const mp_real &eps, double gamma) {
  int n = x.size();
  int history_length = 8;
  int print_interval = 25;
  int check_interval = 100;
  mp_real t;
  mp_real teps = eps * 1.0e20;
  mp_real max_bound = 0.0;

  matrix<mp_real> y(n);
  matrix<mp_real> h(n, n-1);
  matrix<mp_real> b(n, n);
  matrix<mp_real> history(n, history_length);
  bool use_only_one_pair = false;

  int result = RESULT_CONTINUE;
  init_pslq(x, y, h, b);
  history.zero();
  if (debug_level >= 3) {
    y.print("Initial y:");
    b.print("Initial B:");
    h.print("Initial H:");
  }

  pslq_iter = 0;
  while (result == RESULT_CONTINUE) {
    pslq_iter++;
    if (debug_level >= 2) {
      if (debug_level >= 3 || pslq_iter % print_interval == 0)
        cout << "Iteration " << std::setw(5) << pslq_iter << endl;
    }

    result = iterate_mpslq(gamma, y, h, b, eps, teps, use_only_one_pair);
    use_only_one_pair = false;
    
    if (debug_level >= 3) {
      y.print("Updated y: ");
      b.print("Updated B: ");
      h.print("Updated H: ");
    }

    if (result == RESULT_CONTINUE) {
      /* Check the y vector with those of recent iterations. */
      int check_result = check_history(y, history, teps);
      if (check_result == RESULT_DUPLICATE)
        use_only_one_pair = true;
      insert_history(y, history);
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

  if (result == RESULT_RELATION_FOUND) {
    int jm = -1;
    mp_real u, v;

    /* Output final norm bound. */
    jm = matrix_minmax(y, u, v);
    bound_pslq(h, t);
    cout << "Relation detected at iteration " << pslq_iter << endl;
    cout << "  min(y) = " << u << endl;
    cout << "  max(y) = " << v << endl;
    cout << "  bound  = " << t << endl;

    /* Relation found.  Select the relation with smallest norm. */
    if (jm < 0) {
      cerr <<  "ERROR: Invalid index." << endl;
      exit(-1);
    }

    for (int i = 0; i < n; i++) {
      rel(i) = b(jm, i);
    }

  }

  return result;
}

