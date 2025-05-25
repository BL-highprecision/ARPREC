#include <cmath>
#include <cfloat>

#include <arprec/mp_real.h>

#include "pslq1.h"
#include "pslq2.h"


void pslq_dot(int n, int inca, const mp_real *a, 
              int incb, const double *b, mp_real &c) {
  mp_real::mpdotd(n, inca, a, incb, b, c);
}

void pslq_dot(int n, int inca, const mp_real *a, 
              int incb, const mp_real *b, mp_real &c) {
  c = 0.0;
  int ai = 0, bi = 0;
  for (int i = 0; i < n; i++, ai += inca, bi += incb) {
    c += a[ai] * b[bi];
  }
}

