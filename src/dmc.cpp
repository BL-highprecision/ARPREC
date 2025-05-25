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

void mp_real::mpdmc(double a, int n, mp_real& b, int prec_words)
{
  /**
   * This routine converts the DPE number (A, N) to MP form in B.  All bits of
   * A are recovered in B.  However, note for example that if A = 0.1D0 and N
   * is 0, then B will NOT be the multiprecision equivalent of 1/10.  Debug
   * output starts with debug_level = 9.
   *
   */
  int i, k, n1, n2;
  double aa;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(b);
    return;
  }
  
  if (debug_level >= 8) {
    cerr << " MPDMC I: a = " << a << endl;
    cerr << "n = " << n << endl;
    cerr << "mpdmc[1] prec_words " << prec_words << endl;
  }
  
  //  Check for zero.
  if (a == 0.0) {
    zero(b);
    if (debug_level >= 9) print_mpreal("MPDMC O ", b);
    return;
  }
  
  if(n) {
    n1 = n / mpnbt;      // n = mpnbt*n1+n2
    n2 = n - mpnbt * n1;
    aa = ldexp(fabs(a), n2);
  } else {
    n1 = n2 = 0;
    aa = fabs(a);
  }

  //  Reduce AA to within 1 and MPBDX.
  if (aa >= mpbdx) {
    for (k = 1; k <= 100; ++k) {
      aa = mprdx * aa;
      if (aa < mpbdx) {
        n1 = n1 + k;
        break;
      }
    }
  } else if (aa < 1.0) {
    for (k = 1; k <= 100; ++k) {
      aa = mpbdx * aa;
      if (aa >= 1.0) {
        n1 = n1 - k;
        break;
      }
    } 
  }
  
  //  Store successive sections of AA into B.
  double d[8];

  d[2] = n1;
  d[3] = FLOOR_POSITIVE(aa);
  aa = mpbdx * (aa - d[3]);
  d[4] = FLOOR_POSITIVE(aa);
  aa = mpbdx * (aa - d[4]);
  d[5] = FLOOR_POSITIVE(aa);
  d[6] = 0.;
  d[7] = 0.;
  
  for (i = 5; i >= FST_M; --i)
    if (d[i] != 0.) break;
    
  aa = i - 2;
  aa = std::min(aa, (b[0])-5);
  b[1] = sign(aa, a);
  for(i=2;i<int(aa)+FST_M;i++)
    b[i] = d[i];
    
  if (debug_level >= 8) print_mpreal("MPDMC O ", b);
}

