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

void mp_real::mproun(mp_real &a) {
  // This performs rounding and truncation of the MP number A.  It is called
  // by MPNORM, and also by other subroutines when the precision level is
  // reduced by one.  It is not intended to be directly called by the user.
  //
  // Maximum space of A used:  MPNW + 5 cells.
  //
  // The parameter AMX is the absolute value of the largest exponent word
  // allowed for MP numbers.
  // const double amx = 33554432.0;//=pow(2.0, 25); // float: 2.e6
  const double amx = 2.1475e9; //=pow(2.0, 31);
  int i, ia, k, na, n4, AllZero, LoopBreak;
  double a2;

  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(a);
    return;
  }

  // Check for initial zeroes.
  //cerr << "0" << endl;
  a2 = a[2]; // exponent
  a[2] = 0.0;
  ia = a[1] >= 0 ? 1 : -1; // sign (1., a(1))
  na = std::min (int (std::abs (a[1])), prec_words);
  na = std::min (na, int(a[0])-5);
  n4 = na + 4; //index of last addressable word.
  //cerr << "1" << endl;
     
  if (a[FST_M] == 0.) {
    // Find the first nonzero word and shift the entire number left.
    // The length of the result is reduced by the length of the shift.

    AllZero = 1;
    for (i = 4; i <= n4; ++i) {
      if (a[i] != 0.0) {
        AllZero = 0;
        break;
      }
    }
       
    if ( AllZero ) {
      zero(a);
      return;
    }

    k = i - FST_M; // number of leading zeros

    // !dir$ ivdep
    for (i = FST_M; i <= n4 - k; ++i)
      a[i] = a[i+k];

    a2 = a2 - k;
    na -= std::max (k - 2, 0);
    if (k == 2) a[na + FST_M] = 0.0; // BUG FIX
  }
  //cerr << "2" << endl;

  // Perform rounding depending on round_dir.

  if (na == prec_words && round_dir >= 1) {
    //printf("G = %g\n", a[na+3]);
    //cerr << "na+3 = " << na+3 << endl;
    if ( (round_dir == 1 && a[na+3] >= 0.5 * mpbdx) || 
         (round_dir == 2 && a[na+3] >= 1) )
      //cerr << "2.1" << endl;
      //cerr << "a = " << a[na + 3] << endl;
      a[na+2] += 1.0;

    // Zero out unused words.
    a[na+3] = a[na+4] = 0.0;

    // Release carries as far as necessary due to rounding.
    LoopBreak = 0;
    for (i = na + 2; i >= FST_M; --i) {
      if (a[i] < mpbdx) {
        LoopBreak = 1; // goto 140
        break;
      }
      a[i] -= mpbdx;
      ++a[i-1];
    }

    // Release of carries due to rounding continued all the way to the start
    // -- i.e. number was entirely 9's.
    if ( !LoopBreak ) {
      a[3] = a[2];
      na = 1;
      a2++;
    }
  }
  //cerr << "3" << endl;

  // 140
  if (a[na+2] == 0.) {
    // At least the last mantissa word is zero.  Find the last nonzero word
    // and adjust the length of the result accordingly.
    AllZero = 1;
    for (i = na + 1; i >= FST_M; --i) {
      if (a[i] != 0.) {
        AllZero = 0;
        break; // goto 160
      }
    }
    if ( AllZero ) {
      zero(a);
      return;
    }

    // 160
    na = i - 2;
    a[na+4] = 0.0;
  }
  //cerr << "4" << endl;

  // Check for overflow and underflow.
  if (a2 < -amx) {
    if (MPKER[68] != 0) {
      cerr << "*** MPROUN: Exponent underflow." << endl;
      error_no = 68;
      if (MPKER[error_no] == 2)  mpabrt();
    }
  } else if (a2 > amx) {
    if (MPKER[69] != 0) {
      cerr << "*** MPROUN: Exponent overflow." << endl;
      error_no = 69;
      if (MPKER[error_no] == 2)  mpabrt();
    }
  }
  //cerr << "5" << endl;

  // Check for zero.

  if (a[FST_M] == 0.) {
    zero(a);
  } else {
    a[1] = ia >= 0 ? na : -na; // sign (na, ia)
    a[2] = a2;
  }
}

