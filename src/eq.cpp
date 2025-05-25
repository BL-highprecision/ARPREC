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

void mp_real::mpeq(const mp_real& a, mp_real& b, int prec_words)
{
  // This routine sets the MP number B equal to the MP number A.  Debug output
  // starts with debug_level = 10.
  //
  // Max DP space for B: MPNW + 3 cells.
  //
  // The fact that only MPNW + 3 cells, and not MPNW + 4 cells, are copied is
  // important in some routines that increase the precision level by one.
  
  int i, ia, na, nb;
  if (error_no != 0) {
    if (error_no == 99) mpabrt();
    zero(b);
    return;
  }
  if (debug_level >= 10) cerr << "MPEQ" << endl;

  ia = int(sign(1.0, a[1]));
  na = std::min(int(std::abs(a[1])), prec_words);
  nb = std::min(na, int(b[0])-FST_M-1);
  if (na == 0) {
    zero(b);
    return;
  }

  b[1] = sign(nb, ia);
  for (i = 2; i < nb + FST_M; ++i) b[i] = a[i];
  nb = std::min(nb + FST_M + 1, int(b[0]));
  for (; i < nb; i++) b[i] = 0.0;
}

