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

#include <sstream>
#include <string>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;
using std::string;

/* convert std::string to double.  Returns true if successful. */
static bool to_double(const string &str, double &v) {
  std::istringstream is(str);
  return !(is >> v).fail();
}

bool mp_real::construct(const string &expn_str, const string &digit_str1, const string &digit_str2) {
  const unsigned int digits_per_step = 6;
  const int factor_per_step = 1000000;
#if 0
  cerr << "expn_str  = '" << expn_str << '\'' << endl;
  cerr << "digit_str = '" << digit_str1;
  if (digit_str2.length() > 0)
    cerr << '.' << digit_str2;
  cerr << '\'' << endl;
#endif
  string str;
  string digits = digit_str1 + digit_str2;
  double e = 0.0;

  if (expn_str.length() > 0) {
    if (!to_double(expn_str, e)) {
      cerr << "mp_real::construct: error reading exponent" << endl;
      return false;
    }
  }

  unsigned int i = 0;
  string::size_type len = digits.length();
  double v;
  bool neg = false;
  int prec_words = static_cast<int>(mpr[0]) - 5;

  if (digits[0] == '-') { i = 1; neg = true; }

  str = digits.substr(i, digits_per_step);
  i += digits_per_step;
  if (!to_double(str, v)) {
    cerr << "mp_real::construct: error reading digits" << endl;
    return false;
  }

  mpr[1] = 1.0;
  mpr[2] = 0.0;
  mpr[3] = v;

  while (i < len) {
    str = digits.substr(i, digits_per_step);

    if (!to_double(str, v)) {
      cerr << "mp_real::construct: error reading digits" << endl;
      return false;
    }

    string::size_type n = str.length();
    if (n < digits_per_step) {
      int f = 1;
      for (unsigned int j = 0; j < n; j++) f *= 10;
      mpmuld(*this, static_cast<double>(f), 0, *this, prec_words);
    } else {
      mpmuld(*this, static_cast<double>(factor_per_step), 0, *this, prec_words);
    }
    mpadd(*this, mp_real(v, 8), *this, prec_words);
    i += 6;
  }

  if (neg) mpr[1] = -mpr[1];
  mp_real t(10.0, 6), s(0.0, prec_words + 6);
  mpnpwr(t, static_cast<int>(e) - static_cast<int>(digit_str2.length()), s, prec_words);
  mpmul(*this, s, *this, prec_words);
  mproun(*this);

  return true;
}

