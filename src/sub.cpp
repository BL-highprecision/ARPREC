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

void mp_real::mpsub(const mp_real &a, const mp_real &b, mp_real& c,
			   int prec_words)
{
  // This routine subtracts MP numbers A and B to yield the MP difference C,
  // by negating B and adding.  Debug output starts with debug_level = 9.
  //
  // Max SP space for C: MPNW + 5 cells.

  int i, BreakLoop;
  double b1;

  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(c);
    return;
  }
  
  if (debug_level >= 9) cerr << " MPSUB" << endl;

  // Check if A = B.  This is necessary because A and B might be same array,
  // in which case negating B below won't work.

  // check if A == B points to the same object 
  if(&a == &b) {
    zero(c);
    if(debug_level >= 9) print_mpreal("MPSUB O ", c);
    return;
  }
  
  // check if their exponent and mantissas are the same
  if (a[1] == b[1]) {
    BreakLoop = 0;
    for (i = 2; i < int(std::abs(a[1])) + FST_M; ++i) {
      if (a[i] != b[i]) {
        BreakLoop = 1;
        break;
      }
    }
    if (!BreakLoop) {
      zero(c);
      if(debug_level >= 9) print_mpreal("MPSUB O ", c);
      return;
    }
  }
  
  // Save the sign of B, and then negate B.
  b1 = b[1];
  double *temp; // use temp to keep const modifier
  temp = b.mpr;
  temp[1] = -b1;

  // Perform addition and restore the sign of B.
  mpadd(a, b, c, prec_words);
  
  // When restoring the sign of b, we must make sure that
  // b and c were not the same object.  if they were,
  // then b was overwriten, and c already contains the correct
  // result.
  if(&b != &c)
     temp[1] = b1;

  return;
}

