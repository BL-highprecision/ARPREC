#include <iostream>
#include <iomanip>
#include <cmath>

#include <arprec/mp_real.h>
#include <arprec/mp_int.h>

#include "pslq2.h"

using std::cerr;
using std::cout;
using std::endl;

int pslq2(const matrix<mp_real> &x, matrix<mp_real> &rel, 
          const mp_real &eps, double gamma) {
  int n = x.size();
  int print_interval = 250;
  int check_interval = 10;
  int save_interval  = 10;
  mp_real t;
  mp_real teps = eps * 1.0e20;
  double max_bound = 0.0;

  matrix<mp_real> y(n);
  matrix<mp_real> h(n, n-1);
  matrix<mp_real> b(n, n);

  matrix<double> dy(n);
  matrix<double> da(n, n);
  matrix<double> db(n, n);
  matrix<double> dh(n, n-1);

  matrix<double> dsy(n);
  matrix<double> dsa(n, n);
  matrix<double> dsb(n, n);
  matrix<double> dsh(n, n-1);

  int result;

  std::ios_base::fmtflags fmt = cout.flags();

  cout << std::scientific << std::setprecision(20);

  init_pslq(x, y, h, b);
  result = reduce_pslq(h, y, b, eps);

  if (debug_level >= 3) {
    y.print("Initial y:");
    b.print("Initial B:");
    h.print("Initial H:");
  }

  pslq_iter = 0;
  while (result == RESULT_CONTINUE) {

    result = init_pslq2(dy, dh, da, db, y, h, 1.0e-10, 6);
    if (debug_level >= 3) {
      dy.print("Initial dy:");
      da.print("Initial da:");
      db.print("Initial db:");
      dh.print("Initial dh:");
    }
    if (result == RESULT_CONTINUE) {

      int iter_saved = pslq_iter;
      lq_decomp(n, n-1, dh);
      save_pslq(da, db, dh, dy, dsa, dsb, dsh, dsy);

      if (debug_level >= 3) {
        dh.print("after LQ factorization: dh");
      }

      while (result == RESULT_CONTINUE) {

        if (pslq_iter - iter_saved >= save_interval) {
          if (debug_level >= 3)
            cout << "Iteration " << std::setw(6) << pslq_iter << ": DP save" << endl;
          save_pslq(da, db, dh, dy, dsa, dsb, dsh, dsy);
          iter_saved = pslq_iter;
        }

        pslq_iter++;
        if (debug_level >= 2) {
          if (debug_level >= 3 || pslq_iter % print_interval == 0)
            cout << "Iteration " << std::setw(6) << pslq_iter << endl;
        }

        result = iterate_pslq2(gamma, da, db, dh, dy);

        if (debug_level >= 3) {
          dy.print("Updated dy: ");
          db.print("Updated dB: ");
          da.print("Updated dA: ");
          dh.print("Updated dH: ");
        }

      }

      if (result == RESULT_VERY_LARGE_VALUE) {
        /* Double precision iteration resulted in loss of precision.
           Restore to previous data.                                    */
        if (debug_level >= 2) {
          cout << "Iteration " << std::setw(6) << pslq_iter << ": DP iteration aborted" << endl;
          cout << "  Reverting to iteration " << iter_saved << endl;
        }
        pslq_iter = iter_saved;
        save_pslq(dsa, dsb, dsh, dsy, da, db, dh, dy);
      }
    }

    if (result == RESULT_LARGE_VALUE       ||
        result == RESULT_VERY_LARGE_VALUE  ||
        result == RESULT_RELATION_FOUND) {

      /* Update the MP arrays. */
      int update_result = update_pslq(da, db, b, h, y, eps, teps);
      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(6) << pslq_iter << ": MP update" << endl;
      }
      if (debug_level >= 3) {
        y.print("Updated y: ");
        b.print("Updated B: ");
        h.print("Updated H: ");
      }

      /* Compute norm bound. */
      lq_decomp(n, n-1, dh);
      double bnd;
      bound_pslq(dh, bnd);
      if (bnd > max_bound)
        max_bound = bnd;

      if (bnd < 1.0e-300) {
        cerr << "ERROR: bound too small" << endl;
        exit(-1);
      }

      if (debug_level >= 2) {
        mp_real u, v;
        matrix_minmax(y, u, v);

        cout << "Iteration " << std::setw(6) << pslq_iter << endl;
        cout << "  norm bound = " << bnd << endl;
        cout << "  max  bound = " << max_bound << endl;
        cout << "  min(y)     = " << u << endl;
        cout << "  max(y)     = " << v << endl;;
      }

      if (update_result == RESULT_CONTINUE &&
          result != RESULT_VERY_LARGE_VALUE)
        result = RESULT_CONTINUE;
      else
        result = update_result;
    }

    if (result == RESULT_RANGE_TOO_LARGE ||
        result == RESULT_VERY_LARGE_VALUE) {
      /* Perform MP iterations until we can start 
         double-precision iterations. */
      lq_decomp(n, n-1, h);

      int iter_checked = pslq_iter;
      do {
        if (pslq_iter - iter_checked >= check_interval) {
          int check_result = check_pslq(y, 1.0e-10, 6);
          if (check_result == RESULT_CONTINUE) {
            if (debug_level >= 2)
              cout << "Iteration " << std::setw(6) << pslq_iter << ": Return to DP iterations" << endl;
            break;
          }
          iter_checked = pslq_iter;
        }

        pslq_iter++;
        result = iterate_pslq(gamma, y, h, b, eps, teps);
      } while (result == RESULT_CONTINUE);
    }

  } /* outer while */

  int jm = -1;
  if (result == RESULT_RELATION_FOUND) {
    mp_real u, v;

    /* Output final norm bound. */
    /*Relation found.  Select the relation with smallest y and compute norm.*/
    jm = matrix_minmax(y, u, v);
    bound_pslq(h, t);
    cout << "Relation detected at iteration " << pslq_iter << endl;
    cout << "  min(y) = " << u << endl;
    cout << "  max(y) = " << v << endl;
    cout << "  bound  = " << t << endl;

    if (jm < 0) {
      cerr << "ERROR: Invalid index." << endl;
      exit(-1);
    }

    for (int i = 0; i < n; i++) {
      rel(i) = b(i, jm);
    }

  }

  cout.flags(fmt);
  return result;
}

