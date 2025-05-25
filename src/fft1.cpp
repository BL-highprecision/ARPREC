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

void mp_real::mpfft1(int is, int m, int m1, int m2, 
		     double *x, double *y) 
{
  /**
   * 
   *
   * This routine performs a complex-to-complex FFT. IS is the sign of
   * the transform (1 or -1). N = 2^m is the size of the transform.  N1 =
   * 2^M1 and N2 = 2^M2, where M1 and M2 are defined below.  Xi s the input and
   * output array, and Y is a scratch array.  X must have at least N, and Y 
   * at least N + N1 *MPNSP1 double-double complex cells.  The arrays MPUU1
   * and MPUU2 must have beeninitialized by calling MPINIX.  This routine is
   * not intendedto be called directly by the user.
   *
   * This employs the two pass variant of the "four-step" FFT.  See the article
   * by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35. 
   *
   */

  //int n = 1 << m; // 2 ^ m
  //m1 = (m+1)/2;
  int n1 = 1 << m1;
  //m2 = m - m1;
  int n2 = 1 << m2;
  int nr1 = std::min(n1, mpnrow),
    nr2 = std::min(n2, mpnrow);
  int ku = int(mpuu2[2*(m-1)]);
  int i, j, k, iu;
  int n_fft;

  //z1(mpnrow+mpnsp1, n1), x(n1, n2), y(n2+mpnsp1, n1), 
  // z2(mpnrow+mpnsp1, n1)
  // are complex arrays
  int z1_lda =2*(mpnrow)+mpnsp1;
  int x_lda = n1*2;
  int y_lda = 2*(n2)+mpnsp1;
  int z2_lda = z1_lda;
  int z3_lda = z1_lda;
  int z1_col_start, y_col_start, x_col_start, mpuu2_start,
    z2_col_start;
  int z3_col_start;
  double *z1 = (new double[z1_lda*n1]);
  double *z2 = (new double[z2_lda*n1]);
  double *z3;
  
  for(i=0;i<n1;i += nr1) {
    n_fft = std::min(nr1, n1 - i);
    // Copy Nr1 rows of X (treated as a N1 x N2 complex array) into Z1
    for(j=0;j<n2;j++) {
      z1_col_start = z1_lda * j;
      x_col_start = x_lda * j;
      for(k=0;k<n_fft*2;k+=2) {
	z1[z1_col_start+k] = x[x_col_start + (i*2)+k];
	z1[z1_col_start+k+1] = x[x_col_start+(i*2)+k+1];
      }
    }
    
    // Perform roughly NR1 FFTs, each of length N2.

    mpfft2(is, n_fft, m2, n2, z1, z2, z3);
  
    // Multiply the resulting NR1 x N2 complex block by roots of unity
    // and store transposed into the appropreate section of Y.
   
    iu = i+ ku;// -n1 -1;
    iu = iu*2;
    if(is == 1) {
      if(i == 0) {
	for(k=0;k<n_fft;k++) {
	  y[k*y_lda] = z3[2*k];
	  y[k*y_lda+1] = z3[2*k+1];
	}
	for(j=1;j<n2;j++) {
	  mpuu2_start = iu + j*n1*2;
	  z3_col_start = j * z3_lda;
	  y[i*y_lda + j*2] = z3[z3_col_start];
	  y[i*y_lda + j*2+1] = z3[z3_col_start+1];
	  for(k=1;k<n_fft;k++) {
	    //y(j, i+k) = mpuu2(iu+k+j*n1) * z3(k, j)
	    y[(i+k)*y_lda + (j*2)] = 
	      mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k] -
	      mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k + 1];
	    y[(i+k)*y_lda + (j*2) + 1] = 
	      mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k + 1] +
	      mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k];
	  }
	}
      } else {
	for(j=0;j<n2;j++) {
	  mpuu2_start = iu + j*n1*2;
	  z3_col_start = j * z3_lda;
	  for(k=0;k<n_fft;k++) {
	    //y(j, i+k) = mpuu2(iu+k+j*n1) * z3(k, j)
	    y[(i+k)*y_lda + (j*2)] = 
	      mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k] -
	      mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k + 1];
	    y[(i+k)*y_lda + (j*2) + 1] = 
	      mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k + 1] +
	      mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k];
	  }
	}
      }
    } else {
      //inverse transform, cojugate mpuu2
      j = 0;
      if(i == 0) {
	for(k=0;k<n_fft;k++) {
	  y[k*y_lda] = z3[2*k];
	  y[k*y_lda+1] = z3[2*k+1];
	}
	j++;
      }
      for(/*j=j*/;j<n2;j++) {
	mpuu2_start = iu + j*n1*2;
	z3_col_start = j * z3_lda;
	for(k=0;k<n_fft;k++) {
	  //y(j, i+k) = conj(mpuu2(iu+k+j*n1)) * z3(k, j)
	  y[(i+k)*y_lda + (j*2)] = 
	    mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k] +
	    mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k + 1];
	  y[(i+k)*y_lda + (j*2) + 1] = 
	    mpuu2[mpuu2_start+2*k] * z3[z3_col_start + 2*k + 1] -
	    mpuu2[mpuu2_start+2*k+1] * z3[z3_col_start + 2*k];
	}
      }
    }
  }
  
  for(i=0;i<n2;i+= nr2) {
    // Copy NR2 rows of the Y array into Z2.
    
    n_fft = std::min(nr2, n2 - i);

    for(j=0;j<n1;j++) {
      z2_col_start= j * z2_lda;
      y_col_start = j * y_lda;
      for(k=0;k<n_fft*2;k+=2) {
	z2[z2_col_start + k] = y[i*2+k+y_col_start];
	z2[z2_col_start + k + 1] = y[i*2+k+y_col_start + 1];
      }
    }

    // Perform NR2 FFTs, each of length N1.
    mpfft2(is, n_fft, m1, n1, z2, z1, z3);
    
    // Copy NR2 x N1 complex block back into X array. It's a little more
    // complicated if M is odd.
    
    if(!(m & 0x1)) { // if m is even,...
      for(j=0;j<n1;j++) {
	x_col_start = j * x_lda;
	z3_col_start = j * z3_lda;
	for(k=0;k<2*n_fft;k+=2) {
	  x[x_col_start+(i*2)+k] = z3[z3_col_start + k];
	  x[x_col_start+(i*2)+k+1] = z3[z3_col_start + k + 1];
	}
      }
    } else {
      int j2;
      for(j=0;j<n1/2;j++) {
	j2 = 2 * j;
	x_col_start = j * x_lda;
	z3_col_start = j2 * z3_lda;
	for(k=0;k<n_fft*2;k+=2) {
	  x[x_col_start+(i*2)+k] = z3[z3_col_start + k];
	  x[x_col_start+(i*2)+k + 1] = z3[z3_col_start + k + 1];

	  x[x_col_start+(i*2)+k + 2*n2] = z3[z3_col_start + k + z3_lda];
	  x[x_col_start+(i*2)+k + 2*n2 + 1] =z3[z3_col_start + k + z3_lda + 1];
	}
      }
    }
  }

  delete [] z1;
  delete [] z2;
  return;
}

