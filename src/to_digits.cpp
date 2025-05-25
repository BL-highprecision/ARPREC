/*
 * src/to_digits.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2006
 *
 */
#include <cassert>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

int mp_real::to_digits(char *digits, int &expn, int n) const {
  /* Converts the mp_real this into an array of at most n digits, 
   * along with exponent.  It will append a null character at the end.
   * There must be at least n+1 slots in the character array digits.
   * It returns the length string actually written. */
  int len = 0;
  int nn, j;
  const double al2 = 0.301029995663981195; // log10(2)
  const int digits_per_step = 6;
  const int factor_per_step = 1000000;
  double aa, t1;
  int prec_words = static_cast<int>(n / log10(mpbdx)) + 3;
  int n6 = prec_words + 8;
  mp_real f(10.0, 9), sk0(0.0, n6), sk1(0.0, n6);
  
  int na = std::min (static_cast<int>(std::abs(mpr[1])), prec_words);

  if (na == 0) {
    // zero-length mantissa
    digits[len++] = '0';
    digits[len] = 0;
    expn = 0;
    return len;
  }

  // get approximate exponent
  aa = mpr[3];
  if (na >= 2) aa = aa + mprdx * mpr[4];
  if (na >= 3) aa = aa + mprx2 * mpr[5];
  if (na >= 4) aa = aa + mprdx * mprx2 * mpr[6];
  t1 = al2 * mpnbt * mpr[2] + log10 (aa);
  if (t1 >= 0.0) 
    expn = static_cast<int>(t1);
  else 
    expn = static_cast<int>(t1 - 1.0);

  if(expn >= 0) {
    mpnpwr(f, expn - digits_per_step + 1, sk0, prec_words);
    mpdiv(*this, sk0, sk1, prec_words);
  } else {
    mpnpwr(f, -expn + digits_per_step - 1, sk0, prec_words);
    mpmul(*this, sk0, sk1, prec_words);
  }

  // if we didn't quite get it exactly right, multiply or divide by 10 to fix.
  for (;;) {
    assert(sk1[2] == 0.0);
    if (sk1[3] < factor_per_step / 10.0) {
      --expn;
      mpmuld(sk1, 10.0, 0, sk0, prec_words);
      mpeq(sk0, sk1, prec_words);
    } else if (sk1[3] >= factor_per_step) {
      ++expn;
      mpdivd(sk1, 10.0, 0, sk0, prec_words);
      mpeq(sk0, sk1, prec_words);
    } else {
      break;
    }
  }
  sk1[1] = std::abs(sk1[1]);

  // extract digits_per_step digits at a time
  int mpnw_change = digits_per_step * (static_cast<int>(mpnbt * log10(2.0) / digits_per_step) + 2);
  int factor = 0;
  while (len < n) {
    if (sk1[2] == 0.0) {
      assert(sk1[3] < factor_per_step);
      nn = static_cast<int>(sk1[3]);
      
      factor = factor_per_step / 10;
      int k = nn;
      int d;
      while (len < n && factor > 0) {
        d = k / factor;
        digits[len++] = static_cast<char>(d) + '0';
        k -= d * factor;
        factor /= 10;
      }

      f[1] = 1.0; 
      f[3] = nn - k;
      mpsub(sk1, f, sk0, prec_words); 
      if (sk0[1] == 0) break;
      if (len < n) 
        mpmuld(sk0, static_cast<double>(factor_per_step), 0, sk1, prec_words);
    } else {
      factor = 0;
      for (int i = 0; i < digits_per_step && len < n; i++, len++) {
        digits[len] = '0';
      }
      if (len < n)
        mpmuld(sk1, static_cast<double>(factor_per_step), 0, sk1, prec_words);
    }
    if (len % mpnw_change == 0) prec_words--;
  }

  // check if last digits needs to be rounded up
  if (sk0[1] > 0.0 && sk0[2] >= -1.0) {
    aa = sk0[3];
    if (sk0[1] > 1.0) aa += sk0[4] * mprdx;
    if (sk0[2] == -1.0) aa *= mprdx;
    if (factor > 0) 
      aa /= (10 * factor);
    if (aa >= 0.5) {
      int k = len-1;
      digits[k]++;
      while (k >= 1 && digits[k] > '9') {
        digits[k] = '0';
        digits[--k]++;
      }
      if (digits[0] > '9') {
        digits[0] = '1';
        expn++;
      }
    }
  }

  // trim any trailing zeros
  j = len-1;
  while (j > 0 && digits[j] == '0') j--;
  len = j+1;

  digits[len] = 0;
  return len;
}

