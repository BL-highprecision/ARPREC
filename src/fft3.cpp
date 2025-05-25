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

void mp_real::mpfft3_radix2(int is, int m, int ns, int junk1, int junk2,
		     double *x, double *y) 
{
  /**
   * Computes the first iteration of the radix two Stockham
   *  fft.. used only by mpfft2.
   * This routine is not intended to be called by the user
   *
   * arguments is, junk1, and junk2 are ignored. They are included
   * so that the argument list matches that of mpfft3.
   */
  int i, j, n_half, lda = 2*mpnrow +mpnsp1;
  int pos = 0;

  //n = 1 << m;
  n_half = 1 << (m-1);
  int other_half = lda * n_half;
  
  for(j=0;j<n_half;j++) {
    
    for(i=0;i<ns;i++) {
      y[pos+i*2] = x[pos+i*2] + x[pos+other_half+i*2];
      y[pos+i*2+1] = x[pos+i*2+1] + x[pos+other_half+i*2+1];
      
      y[pos+i*2+other_half] = x[pos+i*2] - x[pos+other_half+i*2];
      y[pos+i*2+other_half+1] = x[pos+i*2+1] - x[pos+other_half+i*2+1];
    }
    pos += lda;
  }
  return;
}

void mp_real::mpfft3(int is, int l, int ns, int m, int n, 
		     double *x, double *y) 
{
  /** 
   * This performs the L-th iteration of 
   * the Stockham FFT on the NS vectors in X.  Y is a scratch array, and UU1
   * is the root of unity array.  X, Y, and UU1 are double complex.  This
   * routine is not intended to be called directly by the user.
   */

  // dimensions : x[mpnrow+mpnsp1, n], y[mpnrow+mpnsp1, n].
  int lda = 2 * mpnrow + mpnsp1 ;
  int x_start_1, y_start_1, x_start, y_start;
  int next_quarter, next_quarter_y;
  double u1[2], u2[2], u3[2];
  register double z1[2], z2[2], z3[2], z4[2], z5[2];
  
  // set initial parameters.
  int L = 1 << l;
  int r = n >> l; // = n / L
  int Ls = L /4;
  //int rs = r * 4;
  int j, k, i;
  int rlda = r*lda;
  int rslda = rlda * 4; // = rs * lda
  int ku = 1 << l, ku2 = ku/2;
  next_quarter = rlda;
  next_quarter_y = rlda << (l-2); // = Ls*r*lda;
  if(is == 1) {
    x_start_1 = y_start_1 = 0;
    for(k=0;k<r;k++) {
      for(i=0;i<ns;i++) {
	z1[0] = x[x_start_1 + i*2];
	z1[1] = x[x_start_1 + 1 + i*2];
	
	z2[0] = x[x_start_1+next_quarter + i*2];
	z2[1] = x[x_start_1+next_quarter + 1 + i*2];
	
	z5[0] = z1[0] - x[x_start_1+2*next_quarter + i*2];
	z5[1] = z1[1] - x[x_start_1+2*next_quarter + i*2 + 1];
	
	z1[0] += x[x_start_1+2*next_quarter + i*2];
	z1[1] += x[x_start_1+2*next_quarter + i*2 + 1];
	
	z3[0] = z2[0] - x[x_start_1+3*next_quarter + i*2]; 
	z3[1] = z2[1] - x[x_start_1+3*next_quarter + i*2 + 1];
	
	z2[0] += x[x_start_1+3*next_quarter + i*2];
	z2[1] += x[x_start_1+3*next_quarter + i*2 + 1]; 
	
	y[y_start_1 + i*2] = z1[0] + z2[0];
	y[y_start_1+1 + i*2] = z1[1] + z2[1];
	
	y[y_start_1+next_quarter_y + i*2] = z5[0] - z3[1];
	y[y_start_1+next_quarter_y+1 + i*2] = z5[1] + z3[0];
	y[y_start_1+3*next_quarter_y + i*2] = z5[0] + z3[1];
	y[y_start_1+3*next_quarter_y+1 + i*2] = z5[1] - z3[0];
	
	y[y_start_1+2*next_quarter_y + i*2] = z1[0] - z2[0];
	y[y_start_1+2*next_quarter_y+ 1 + i*2] = z1[1] - z2[1];
	
      }
      x_start_1 += lda;
      y_start_1 += lda;
    } 
    x_start = 0;
    y_start = 0;
    for(j=1;j<Ls;j++) {
      u1[0] = mpuu1[ku+j*2];
      u1[1] = mpuu1[ku+j*2+1];
      u2[0] = mpuu1[ku2+j*2];
      u2[1] = mpuu1[ku2+j*2+1];
      u3[0] = mpuu3[ku+j*2];
      u3[1] = mpuu3[ku+j*2+1];
      x_start += rslda; // = j*rslda;
      y_start += rlda;// =  j*rlda;
    
      x_start_1 = x_start;
      y_start_1 = y_start;

      for(k=0;k<r;k++) {
	for(i=0;i<ns;i++) {
	  z1[0] = x[x_start_1 + i*2];
	  z1[1] = x[x_start_1 + 1 + i*2];
	  z2[0] = u1[0] * x[x_start_1+next_quarter + i*2] - 
	    u1[1] * x[x_start_1+next_quarter +1 + i*2];
	  z2[1] = u1[0] * x[x_start_1+next_quarter + 1 + i*2] +
	    u1[1] * x[x_start_1+next_quarter + i*2];
	  
	  z3[0] = u2[0] * x[x_start_1+2*next_quarter + i*2] - 
	    u2[1] * x[x_start_1+2*next_quarter +1 + i*2];
	  z3[1] = u2[0] * x[x_start_1+2*next_quarter + 1 + i*2] +
	    u2[1] * x[x_start_1+2*next_quarter + i*2];
	  
	  z4[0] = u3[0] * x[x_start_1+3*next_quarter + i*2] - 
	    u3[1] * x[x_start_1+3*next_quarter +1 + i*2];
	  z4[1] = u3[0] * x[x_start_1+3*next_quarter + 1 + i*2] +
	    u3[1] * x[x_start_1+3*next_quarter + i*2];
	  
	  z5[0] = z1[0] - z3[0];
	  z5[1] = z1[1] - z3[1];
	
	  z1[0] += z3[0];
	  z1[1] += z3[1];
	
	  z3[0] = z2[0] - z4[0];
	  z3[1] = z2[1] - z4[1];

	  z2[0] += z4[0];
	  z2[1] += z4[1];
	
	  y[y_start_1 + i*2] = z1[0] + z2[0];
	  y[y_start_1+1 + i*2] = z1[1] + z2[1];

	  y[y_start_1+next_quarter_y + i*2] = z5[0] - z3[1];
	  y[y_start_1+next_quarter_y+1 + i*2] = z5[1] + z3[0];
	  y[y_start_1+3*next_quarter_y + i*2] = z5[0] + z3[1];
	  y[y_start_1+3*next_quarter_y+1 + i*2] = z5[1] - z3[0];
	  
	  y[y_start_1+2*next_quarter_y + i*2] = z1[0] - z2[0];
	  y[y_start_1+2*next_quarter_y+ 1 + i*2] = z1[1] - z2[1];
	
	}
	x_start_1 += lda;
	y_start_1 += lda;
      } 
    }
  } else {
    x_start_1 = y_start_1 = 0;
    for(k=0;k<r;k++) {
      for(i=0;i<ns;i++) {
	z1[0] = x[x_start_1 + i*2];
	z1[1] = x[x_start_1 + 1 + i*2];
	
	z2[0] = x[x_start_1+next_quarter + i*2];
	z2[1] = x[x_start_1+next_quarter + 1 + i*2];
	
	z5[0] = z1[0] - x[x_start_1+2*next_quarter + i*2];
	z5[1] = z1[1] - x[x_start_1+2*next_quarter + i*2 + 1];
	
	z1[0] += x[x_start_1+2*next_quarter + i*2];
	z1[1] += x[x_start_1+2*next_quarter + i*2 + 1];
	
	z3[0] = z2[0] - x[x_start_1+3*next_quarter + i*2]; 
	z3[1] = z2[1] - x[x_start_1+3*next_quarter + i*2 + 1];
	
	z2[0] += x[x_start_1+3*next_quarter + i*2];
	z2[1] += x[x_start_1+3*next_quarter + i*2 + 1]; 
	
	y[y_start_1 + i*2] = z1[0] + z2[0];
	y[y_start_1+1 + i*2] = z1[1] + z2[1];
	
	y[y_start_1+next_quarter_y + i*2] = z5[0] + z3[1];
	y[y_start_1+next_quarter_y+1 + i*2] = z5[1] - z3[0];
	y[y_start_1+3*next_quarter_y + i*2] = z5[0] - z3[1];
	y[y_start_1+3*next_quarter_y+1 + i*2] = z5[1] + z3[0];
	
	y[y_start_1+2*next_quarter_y + i*2] = z1[0] - z2[0];
	y[y_start_1+2*next_quarter_y+ 1 + i*2] = z1[1] - z2[1];
	
      }
      x_start_1 += lda;
      y_start_1 += lda;
    } 
    x_start = 0;
    y_start = 0;
    for(j=1;j<Ls;j++) {
      u1[0] = mpuu1[ku+j*2];
      u1[1] = mpuu1[ku+j*2+1];
      u2[0] = mpuu1[ku2+j*2];
      u2[1] = mpuu1[ku2+j*2+1];
      u3[0] = mpuu3[ku+j*2];
      u3[1] = mpuu3[ku+j*2+1];
      x_start += rslda; // = j*rslda;
      y_start += rlda;  // = j*rlda;
      x_start_1 = x_start;
      y_start_1 = y_start;


      for(k=0;k<r;k++) {
	for(i=0;i<ns;i++) {
	  z1[0] = x[x_start_1 + i*2];
	  z1[1] = x[x_start_1 + 1 + i*2];
	
	  //conjugate roots of unity and multiply
	  z2[0] = u1[0] * x[x_start_1+next_quarter + i*2] + 
	    u1[1] * x[x_start_1+next_quarter +1 + i*2];
	  z2[1] = u1[0] * x[x_start_1+next_quarter + 1 + i*2] -
	    u1[1] * x[x_start_1+next_quarter + i*2];
	  
	  z3[0] = u2[0] * x[x_start_1+2*next_quarter + i*2] + 
	      u2[1] * x[x_start_1+2*next_quarter +1 + i*2];
	  z3[1] = u2[0] * x[x_start_1+2*next_quarter + 1 + i*2] -
	    u2[1] * x[x_start_1+2*next_quarter + i*2];
	  
	  z4[0] = u3[0] * x[x_start_1+3*next_quarter + i*2] + 
	    u3[1] * x[x_start_1+3*next_quarter +1 + i*2];
	  z4[1] = u3[0] * x[x_start_1+3*next_quarter + 1 + i*2] -
	    u3[1] * x[x_start_1+3*next_quarter + i*2];

	  z5[0] = z1[0] - z3[0];
	  z5[1] = z1[1] - z3[1];
	
	  z1[0] += z3[0];
	  z1[1] += z3[1];
	
	  z3[0] = z2[0] - z4[0];
	  z3[1] = z2[1] - z4[1];

	  z2[0] += z4[0];
	  z2[1] += z4[1];
	
	  y[y_start_1 + i*2] = z1[0] + z2[0];
	  y[y_start_1+1 + i*2] = z1[1] + z2[1];

	  y[y_start_1+next_quarter_y + i*2] = z5[0] + z3[1];
	  y[y_start_1+next_quarter_y+1 + i*2] = z5[1] - z3[0];
	  y[y_start_1+3*next_quarter_y + i*2] = z5[0] - z3[1];
	  y[y_start_1+3*next_quarter_y+1 + i*2] = z5[1] + z3[0];
      
	  y[y_start_1+2*next_quarter_y + i*2] = z1[0] - z2[0];
	  y[y_start_1+2*next_quarter_y+ 1 + i*2] = z1[1] - z2[1];
	
	}
	x_start_1 += lda;
	y_start_1 += lda;
      } //end for
    } //end for
      
  }//end if(is == 1) {} else {}
}

