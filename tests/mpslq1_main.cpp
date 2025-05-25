#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <arprec/mp_real.h>

#include "mpslq1.h"
#include "pslq_main.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) {
  int mode = 0;
  int n;
  int r = 5, s = 5;
  int nr_digits   = 180;
  int n_eps;

  /* Parse command line arguments. */
  parse_command(argc, argv, mode, n, r, s, nr_digits, n_eps);

  n = r * s + 1;
  n_eps = (nr_digits < 700 ? 10 : 20) - nr_digits;

  cout << "nr_digits = " << nr_digits << endl;
  cout << "debug_level = " << debug_level << endl;
  cout << "n = " << n << endl;
  if (debug_level > 0) {
    cout << "r = " << r << "    s = " << s << endl;
    cout << "n_eps = " << n_eps << endl;;
  }

  /* Initialize data */
  mp::mp_init(nr_digits);
  matrix<mp_real> x(n);
  matrix<mp_real> rel(n);
  mp_real eps = pow(mp_real(10.0), n_eps);

  init_data(mode, n, r, s, x, rel);

  if (debug_level > 0) {
    x.print("Initial x:");
  }

  /* Perform Level-1 PSLQ. */
  int result = mpslq1(x, rel, eps);

  /* Output recovered relation. */
  if (result == RESULT_RELATION_FOUND) {
    cout << "Relation found:" << endl;
    cout << std::fixed << std::setprecision(0);
    for (int i = 0; i < n; i++) {
      cout << std::setw(3) << i;
      cout << std::setw(24) << rel(i) << endl;
    }
  } else {
    cout << "Precision exhausted." << endl;
  }

  mp::mp_finalize();
  return 0;
}

