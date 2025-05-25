#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <arprec/mp_real.h>

#include "pslq3.h"

using std::cerr;
using std::cout;
using std::endl;

int pslq3(const matrix<mp_real> &x, matrix<mp_real> &rel, 
          const mp_real &eps, double gamma) {
  int n = x.size();
  int nw  = mp::n_words;
  int nw2 = 10;
  int print_interval = 250;
  int check_interval = 10;
  int save_interval  = 10;
  mp_real t;
  mp_real teps = eps * 1.0e20;
  mp_real max_bound = 0.0;

  matrix<mp_real> y(n);
  matrix<mp_real> h(n, n-1);
  matrix<mp_real> b(n, n);

  SET_PREC(nw2);
  mp_real wteps = pow(mp_real(2.0), 72 - 48 * nw2);
  mp_real weps = wteps * ldexp(1.0, -72);
  mp_real wmx1 = pow(mp_real(2.0), 48 * nw2 - 72);
  mp_real wmx2 = pow(mp_real(2.0), 48 * nw2);

  matrix<mp_real> wy(n);
  matrix<mp_real> wh(n, n-1);
  matrix<mp_real> wa(n, n);
  matrix<mp_real> wb(n, n);
  SET_PREC(nw);

  matrix<double> dy(n);
  matrix<double> da(n, n);
  matrix<double> db(n, n);
  matrix<double> dh(n, n-1);

  matrix<double> dsy(n);
  matrix<double> dsa(n, n);
  matrix<double> dsb(n, n);
  matrix<double> dsh(n, n-1);

  int result;
  init_pslq(x, y, h, b);
  result = reduce_pslq(h, y, b, eps);

  if (debug_level >= 3) {
    y.print("Initial y:");
    b.print("Initial B:");
    h.print("Initial H:");
  }

  pslq_iter = 0;
  if (debug_level >= 2) {
    if (result == RESULT_CONTINUE)
      cout << "Iteration " << std::setw(6) << pslq_iter << ": Starting MP iterations." << endl;
  }
  while (result == RESULT_CONTINUE) {

    result = init_pslq2(wy, wh, wa, wb, y, h, wteps, nw2);
    if (debug_level >= 2) {
      if (result == RESULT_CONTINUE)
        cout << " Iteration " << std::setw(6) << pslq_iter << "Starting MPM iterations." << endl;
    }
    while (result == RESULT_CONTINUE) {

      result = init_pslq2(dy, dh, da, db, wy, wh, 1.0e-10, 5);
      if (result == RESULT_CONTINUE) {

        int iter_saved = pslq_iter;
        lq_decomp(n, n-1, dh);
        save_pslq(da, db, dh, dy, dsa, dsb, dsh, dsy);

        if (debug_level >= 3)
          dh.print("after LQ factorization: dh");

        if (debug_level >= 2)
          cout << "Iteration " << std::setw(6) << pslq_iter << ": Starting DP iterations." << endl;
        while (result == RESULT_CONTINUE) {

          if (pslq_iter - iter_saved >= save_interval) {
            if (debug_level >= 3)
              cout << "Iteration " << std::setw(6) << pslq_iter << ": DP save." << endl;
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

        } /* end while */

        if (result == RESULT_VERY_LARGE_VALUE) {
          /* Double precision iteration resulted in loss of
             precision.  Restore to previously saved data.  */
          if (debug_level >= 2) {
            cout << "Iteration " << std::setw(6) << pslq_iter << ": DP iteration aborted." << endl;
            cout << "                 Reverting to iteration " << iter_saved << endl;
          }
          pslq_iter = iter_saved;
          save_pslq(dsa, dsb, dsh, dsy, da, db, dh, dy);
        }

      } /* end if */

      if (result == RESULT_LARGE_VALUE      ||
          result == RESULT_VERY_LARGE_VALUE ||
          result == RESULT_RELATION_FOUND) {
        /* Update MPM arrays. */
        if (debug_level >= 2) {
          cout << "Iteration " << std::setw(6) << pslq_iter << ": MPM update" << endl;
        }
        int update_result = 
          update_pslq2(da, db, wa, wb, wh, wy, weps, wteps);
        int check_result  = check_size_pslq(wa, wb, wmx1, wmx2);

        mp_real u, v;
        matrix_minmax(wy, u, v);
        if (debug_level >= 2) {
          cout << "Relation detected at iteration " << pslq_iter << endl;
          cout << "  min(wy) = " << u << endl;
          cout << "  max(wy) = " << v << endl;
        }

        if (check_result == RESULT_VERY_LARGE_VALUE) {
          cerr << "ERROR: MPM Array too large." << endl;
          exit(-1);
        }
        if (check_result != RESULT_CONTINUE)
          update_result = check_result;

        if (debug_level >= 3) {
          wy.print("Updated wy: ");
          wa.print("Updated wA: ");
          wb.print("Updated wB: ");
          wh.print("Updated wH: ");
        }

        if (update_result == RESULT_CONTINUE &&
            result != RESULT_VERY_LARGE_VALUE)
          result = RESULT_CONTINUE;
        else
          result = update_result;
      } /* end if */

      SET_PREC(nw2);
      if (result == RESULT_RANGE_TOO_LARGE ||
          result == RESULT_VERY_LARGE_VALUE) {
        if (debug_level >= 2) {
          cout << "Iteration " << std::setw(6) << pslq_iter << ": Performing MPM iterations." << endl;
        }
        /* Perform MPM iterations until we can start 
           double-precision iterations. */
        lq_decomp(n, n-1, wh);
        int iter_checked = pslq_iter;
        do {
          if (pslq_iter - iter_checked >= check_interval) {
            int check_result = check_pslq(wy, 1.0e-10);
            if (check_result == RESULT_CONTINUE) {
              if (debug_level >= 2)
                cout << "Iteration " << std::setw(6) << pslq_iter << ": Return to DP iterations." << endl;
              break;
            }
            iter_checked = pslq_iter;
          }

          pslq_iter++;
          result = iterate_pslq2(gamma, wa, wb, wh, wy, wteps, wmx1, wmx2);
        } while (result == RESULT_CONTINUE);
      }
      SET_PREC(nw);

    } /* end while */

    if (result == RESULT_RANGE_TOO_LARGE ||
        result == RESULT_VERY_LARGE_VALUE) {
      cerr << "MPM range too large." << endl;
      exit(-1);
    }

    if (result == RESULT_LARGE_VALUE ||
        result == RESULT_PRECISION_EXHAUSTED ||
        result == RESULT_RELATION_FOUND) {
      /* Update the MP arrays. */
      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(6) << pslq_iter << ": MP update" << endl;
      }
      int update_result = update_pslq(wa, wb, b, h, y, eps, teps);
      //debug_level = 3;
      mp_real u, v;
      matrix_minmax(y, u, v);
      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(6) << pslq_iter << endl;
        cout << "  min(y) = " << u << endl;
        cout << "  max(y) = " << v << endl;
      }
      if (debug_level >= 3) {
        y.print("Updated y: ");
        b.print("Updated B: ");
        h.print("Updated H: ");
      }

      /* Compute norm bound. */
      SET_PREC(nw2);
      lq_decomp(n, n-1, wh);
      mp_real bnd;
      bound_pslq(wh, bnd);
      if (bnd > max_bound)
        max_bound = bnd;
      SET_PREC(nw);

      if (debug_level >= 2) {
        cout << "Iteration " << std::setw(6) << pslq_iter << endl;
        cout << "  norm bound = " << bnd << endl;
        cout << "   max bound = " << max_bound << endl;
      }

      result = update_result;
    } /* end if */

  } /* end while */

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

  return result;
}

