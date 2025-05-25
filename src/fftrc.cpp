/*
 * src/mprealx.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * This file contains the basic routines for very high levels of precision,
 * > 1000 digits.
 *
 * The main routines among these  mpmulx, which employs FFT for the
 * convolution in multiplication, and mpinix, which initializes certain
 * root of unity arrays for the fft.
 * All other routines in this file are here only to help mpmulx.
 *
 */
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpfftrc(int is, int m, int n, int n_nonzero, double *x, double *y)
{
  /**
   * This performs an N-point real-to-complex FFT, where N == 2^M.  X
   * is the double precision input array, and Y is the double complex output
   * array.  X must be dimensioned with N DP cells, and Y with N/2+N1*NSP1+1
   * cells, where N1 = 2^int(M/2).  This space requirement for Y
   *  is somewhat larger than is used in this function, because
   * MPFFT1, which is called by this function, requires more. IS is the
   * sign of the transform.  Before calling MPFFTRC, the UU1 and UU2 arrays 
   * must initialized by calling mp_initx.  This routine is not intended
   * to be called directly by the user.
   */
  double z1[2], z2[2], a1[2], a2[2], z3[2];
  int mx = int(mpuu1[0]), n1, n2,  n4, k, ku;

  // Check if input parameters are invalid.
  
  if((is != 1 && is != -1) || m<3 || m > mx) {
    if(MPKER[67] != 0) {
      cerr <<"*** MPFFTRC: either the UU arrays have not been initialized," << endl;
      cerr <<"             or else one of the input parameters is invalid." << endl;
      error_no = 67;
      if(MPKER[error_no] == 2) mpabrt();
    }
    return;
  }
  
  n1 = 1 << (m /2);
  n2 = n / 2;
  n4 = n / 4;
  double *dc1 = x;
  
  int m1, m2;


  for(k=n_nonzero;k<n;k++) x[k] = 0.0;
  m1 = ((m-1)+1)/2;
  m2 = (m-1) - m1;

  mpfft1(is, m-1, m1, m2, dc1, y);

  // Reconstruct the FFT of X.
  y[0] = 2.0 * (dc1[0] + dc1[1]);
  y[1] = 0.0;

  if(is == 1) {
    y[n2] = 2.0 * dc1[n2]; y[n2+1] = 2.0 * dc1[n2+1];
  } else {
    y[n2] = 2.0 * dc1[n2]; y[n2+1] = -2.0 * dc1[n2+1];
  }
  y[n] = 2.0 * (dc1[0] - dc1[1]); y[n+1] = 0.0;
  ku = n;

  if(is == 1) {
    for(k=2;k<n2;k+=2) {
      z1[0] = dc1[k]; z1[1] = dc1[k+1];
      z2[0] = dc1[n-k]; z2[1] = -dc1[n+1-k];
      //a1 = z1 + z2
      a1[0] = z1[0] + z2[0]; a1[1] = z1[1] + z2[1];
      
      z3[0] = z1[0] - z2[0]; z3[1] = z1[1] - z2[1];
      //a2 = -i * mpuu1(k+ku) * (z1 - z2)
      a2[0] = (mpuu1[(k+ku)] * z3[1]) + (mpuu1[(k+ku)+1] * z3[0]);
      a2[1] = mpuu1[(k+ku)+1] * z3[1] - mpuu1[(k+ku)] * z3[0];
      
      y[k] = a1[0] + a2[0]; y[k+1] = a1[1] + a2[1];
      //y[n2*2-k] = conj(a1 - a2);
      y[n-k] = a1[0] - a2[0]; y[n-k+1] = a2[1] - a1[1];
    }
  } else {
    for(k=2;k<n2;k+=2) {
      z1[0] = dc1[k]; z1[1] = dc1[k+1];
      z2[0] = dc1[n-k]; z2[1] = dc1[n-k+1];
      //a1 = z1 + z2;
      a1[0] = z1[0] + z2[0]; a1[1] = z1[1] + z2[1];

      //a2 = -i * conj( mpuu1[k+ku])  * (z1 - z2);
      z3[0] = z1[0] - z2[0]; z3[1] = z1[1] - z2[1];
      a2[0] = (mpuu1[k+ku] * z3[1]) - (mpuu1[k+ku+1] * z3[0]);
      a2[1] = - (mpuu1[k+ku] * z3[0]) - (mpuu1[k+ku+1] * z3[1]);
      
      y[k] = a1[0] + a2[0]; y[k+1] = a1[1] + a2[1];
      y[n-k] = a1[0] - a2[0]; y[n-k+1] = a2[1] - a1[1];//conj a1 - a2
    }
  }

  return;    
}

