#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>

mp_real *mp_real::alloc_array(int n) {
  int old_mpnw = prec_words;
  prec_words = 0;
  
  mp_real *m = new mp_real[n];
  int sz = old_mpnw + 5;
  double *raw_data = new double[n * sz];
#if (ARPREC_DEBUG)
  for (int i = 0; i < n*sz; i++) raw_data[i] = mp::_d_nan;
  VALGRIND_MAKE_MEM_UNDEFINED(raw_data, n * sz * sizeof(double));
#endif
  double *p = raw_data;

  for (int i = 0; i < n; i++, p += sz) {
    p[0] = sz;
    p[1] = old_mpnw;
    m[i].mpr = p;
  }

  prec_words = old_mpnw;
  return m;
}

void mp_real::free_array(mp_real *m) {
  double *raw_data = m[0].mpr;
  delete [] m;
  delete [] raw_data;
}

mp_complex *mp_complex::alloc_array(int n) {
  int old_mpnw = prec_words;
  prec_words = 0;

  mp_complex *m = new mp_complex[n];
  int sz = old_mpnw + 5;
  double *raw_data = new double[2 * n * sz];
#if (ARPREC_DEBUG)
  for (int i = 0; i < 2*n*sz; i++) raw_data[i] = mp::_d_nan;
  VALGRIND_MAKE_MEM_UNDEFINED(raw_data, n * sz * sizeof(double));
#endif
  double *p = raw_data;
  for (int i = 0; i < n; i++) {
    p[0] = sz;
    p[1] = old_mpnw;
    m[i].real.mpr = p;
    p += sz;
    p[0] = sz;
    p[1] = old_mpnw;
    m[i].imag.mpr = p;
    p += sz;
  }

  prec_words = old_mpnw;
  return m;
}

void mp_complex::free_array(mp_complex *m) {
  double *raw_data = m[0].real.mpr;
  delete [] m;
  delete [] raw_data;
}

