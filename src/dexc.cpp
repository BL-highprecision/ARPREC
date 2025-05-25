/*
 * src/mpreal4.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Additional routines, mostly those for which speed is less important.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

void mp_real::mpdexc(const char a[], int k, mp_real& b)
{
  /**
   * This routine converts the character*1 string A, which
   * represents a multiprecision number in Fortran style, i.e.
   * '1234567890' or '1.23456789D-21', into standard MP binary format.
   * This routine is not intended to be called directly by the user.
   *
   * This routine has not been thoroughly tested.
   */

  int i;
  int foundExponent = 0;

  for( i = 0; i < k; i++) {
    if (a[i] == 'D' || a[i] == 'E' || a[i] == 'd' ||
      a[i] == 'e') {
      foundExponent = 1;
      break;
    }
  }

  if (!foundExponent) {
    mpinpc (a, k, b);
    return;
  }
	
  char *c = new char[n_digits+101];
  int i1 = i + 1; 
  int k1 = i; 
  int k2 = k - i1; 
  c[0] = '1';
  c[1] = '0';
  c[2] = '^';
  
  for (i = 0; i < k2; ++i) c[i+3] = a[i+i1];
  
  c[k2+3] = 'x';
  
  for (i = 0; i < k1; ++i) c[i+k2+4] = a[i];
  c[i+k2+4]='\0';
	
  mpinpc(c, k1+k2+4, b);
  delete [] c;
  return;
}

