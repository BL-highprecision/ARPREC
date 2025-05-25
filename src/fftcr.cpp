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

void mp_real::mpfftcr(int is, int m, int n, int nsq, 
		      double *x, double *y)
{
  /**
   * This performs an N-point complex-to real FFT, where N = 2^M.
   * X is the double complex input array, and y is the double (real)
   * output array.  The array X is used as a scratch array in MPFFT1, 
   * and so is overwritten.  X must be dimensioned with N/2+N1*NSP1+1 
   * double complex cells, and Y with N double cells, where N = 2^M,
   * and N1 = 2^int(M/2).  this dimension requirement for X is
   * somewhat greater than is used in this function, because MPFFT1, 
   * which is called by this routine, requires more.  IS is the sign of
   * the transform.  Before calling MPFFTCR, the UU1 and UU2 arrays must
   * bt initialized by calling MPINIX.  This routine is not intended 
   * to be called direcly by the user.
   */
  double  a1[2], a2[2], x1[2], x2[2];
  int mx = int(mpuu1[0]), n1, n2, n21, n4, k, ku;

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
  
  n1 = 1 << (m /2 );
  n2 = n / 2;
  n21 = n2 + 1;
  n4 = n / 4;

  // Construct the input to MPFFT1;
  
  y[0] = 0.5 * (x[0] + x[n]);
  y[1] = 0.5 * (x[0] - x[n]);
  
  if(is == 1) {
    y[n2] = x[n2];
    y[n2+1] = -x[n2+1];
  } else {
    y[n2] = x[n2];
    y[n2+1] = x[n2+1];
  }
  ku = n;
  
  if(is == 1) {
    for(k=2;k<n2; k+=2) {
      x1[0] = x[k]; x1[1] = x[k+1];
      
      x2[0] = x[n-k]; x2[1] = -x[n-k+1];
      
      a1[0] = x1[0] + x2[0]; a1[1] = x1[1] + x2[1];
      
      //x1 = x1 - x2;
      x1[0] -= x2[0]; x1[1] -= x2[1];

      // a2 = i * mpuu1[k+ku] * (x1 - x2)
      a2[0] = -(mpuu1[k+ku] * x1[1]) - (mpuu1[k+ku+1] * x1[0]);
      a2[1] = mpuu1[k+ku] * x1[0] - (mpuu1[k+ku+1] * x1[1]);
      
      y[k] = 0.5 * (a1[0] + a2[0]); 
      y[k+1] = 0.5 * (a1[1] + a2[1]);
      
      //y[n-k] = 0.5 * conj(a1-a2)
      y[n-k] = 0.5 * (a1[0] - a2[0]);
      y[n-k+1] = 0.5 * (a2[1] - a1[1]);
    }
  } else {
    for(k=2;k<n2; k+=2) {
      x1[0] = x[k]; x1[1] = x[k+1];
      
      x2[0] = x[n-k]; x2[1] = -x[n-k+1];
      
      a1[0] = x1[0] + x2[0]; a1[1] = x1[1] + x2[1];
      
      //x1 = x1 - x2;
      x1[0] -= x2[0]; x1[1] -= x2[1];

      // a2 = i * conj(mpuu1[k+ku]) * (x1 - x2)
      a2[0] = (mpuu1[k+ku+1] * x1[0]) - (mpuu1[k+ku] * x1[1]);
      a2[1] = (mpuu1[k+ku+1] * x1[1]) + (mpuu1[k+ku] * x1[0]);
      
      y[k] = 0.5 * (a1[0] + a2[0]); 
      y[k+1] = 0.5 * (a1[1] + a2[1]);
      
      //y[n-k] = 0.5 * conj(a1-a2)
      y[n-k] = 0.5 * (a1[0] - a2[0]);
      y[n-k+1] = 0.5 * (a2[1] - a1[1]);
    }
  }
  
  int m1, m2;
  m1 = ((m-1)+1)/2;
  m2 = (m-1) - m1;

  // Perform a normal N/2 point FFT on y.
  mpfft1(is, m-1, m1, m2, y, x);
  
  return;
}


