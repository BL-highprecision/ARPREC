/*
 * src/mpreal.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpmdc(const mp_real&a, double &b, int &n, int prec_words)
{
  /**
   * This procedure takes the mp_real A, and splits it into 
   * a double, b, and a exponent, n. 
   *
   * On exit, the following should be roughly true: 
   *
   *       a ==(roughly) b*2^n
   */

  double aa;
  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    b = 0.0;
    n = 0;
    return;
  }
  if(debug_level >= 9) {
    int no = std::min(int(std::abs(a[1])), debug_words) + 2;
    cerr << "MPMDC I " << no << endl;
  }
  if(a[1] == 0.0) {
    b = 0.0;
    n = 0;
    return;
  }

  int na = int(std::abs(a[1]));
  aa = a[FST_M];
  if(na >= 2) aa += mprdx * a[FST_M+1];
  if(na >= 3) aa += mprx2 * a[FST_M+2];
  if(na >= 4) aa += mprx2 * mprdx * a[FST_M+3];

  n = int(mpnbt * a[2]); 
  b = sign(aa, a[1]);
  
  if(debug_level >= 9) cerr << "MPMDC 0 " << b << ", " << n << endl;
}

