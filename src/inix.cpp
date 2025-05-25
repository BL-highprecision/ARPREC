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

void mp_real::mpinix(int n) 
{
  /**
   * This initializes the root of unity arrays UU1 and UU2, which are 
   * required by the FFT routines that are called by mpmulx.  Before
   * calling any of the advanced MP routines (i.e. those whose names end
   * in X), this routine must be called with N set to the largest precision
   * level prec_words that will be used in the subsequent application.  It is 
   * not necessary for the user to call mpinix if the advanced routines
   * are not called.
   */

  const double cl2 = 1.4426950408889633;
  const double pi = acos(-1.0);
  n *= 2;
  double l_temp = 0.75 * n;
  int m = int(cl2 * log(l_temp) + 1.0 - mprxx);
  int mq = m + 2;
  int nq = 1 << mq;
  if(mpuu1)
    delete [] mpuu1;
  mpuu1 = new double[2 * nq + 2];
  if(mpuu2)
    delete [] mpuu2;
  mpuu2 = new double[2 * (mq + nq + 1)];
  if(mpuu3)
    delete [] mpuu3;
  mpuu3 = new double[2 * nq + 2];
  mpuu1[0] = mq;
  mpuu1[1] = 0.0;
  int ku;
  int ln;
  int i, j;
  double t1, ti, sin_temp, cos_temp;
  
  ku = 2;
  ln = 2;
  mpuu1[2] = 1.0;
  mpuu1[3] = 0.0;
  //mpuu1 holds roots of unity from 0 to pi
  for(j=0;j<mq-1;j++) {
    t1 = pi / double(ln);
    
    for (i = 0; i<ln/2; i+=2) {
      mpuu1[2 * (i+ku)] = mpuu1[i+ku];
      mpuu1[2 * (i+ku) + 1] = mpuu1[i+ku+1];
      ti = double(i+1) * t1;
      sin_temp = sin(ti);
      cos_temp = cos(ti);
      mpuu1[2 * (i+ku+1)] = cos_temp;
      mpuu1[2 * (i+ku+1) + 1] = sin_temp;

      //multiply by i and use again for the roots above pi/2
      mpuu1[2* (i+ku) + ln] = -mpuu1[i+ku+1];
      mpuu1[2* (i+ku) + ln+1] = mpuu1[i+ku]; 

      mpuu1[2* (i+ku) + ln+2] = -sin_temp;
      mpuu1[2* (i+ku) + ln+3] = cos_temp;
      
    }
    
    ku += ln;
    ln *= 2;
  }


  ku = 2;
  ln = 2;
  mpuu3[2] = 1.0;
  mpuu3[3] = 0.0;
  for(j=0;j<mq-1;j++) {
    t1 = pi / double(ln);
    
    for (i = 0; i<ln/* /2*/; i+=2) {
      mpuu3[2 * (i+ku)] = mpuu3[i+ku];
      mpuu3[2 * (i+ku) + 1] = mpuu3[i+ku+1];
      ti = double(3 * (i+1)) * t1;
      sin_temp = sin(ti);
      cos_temp = cos(ti);
      mpuu3[2 * (i+ku+1)] = cos_temp;
      mpuu3[2 * (i+ku+1) + 1] = sin_temp;
    }
    
    ku += ln;
    ln *= 2;
  }
  
  double tpn;
  int nn, nn1, nn2, iu;
  int mm, mm1, mm2, k;

  ku = mq + 1;
  mpuu2[0] = mq;
  
  for(k=1;k<mq-1;k++) {
    mpuu2[2*k] = double(ku-1);
    mpuu2[2*k+1] = double(0.0);
    mm = k + 1;
    nn = 1 << mm;
    mm1 = (mm + 1) / 2;
    mm2 = mm - mm1;
    nn1 = 1 << mm1;
    nn2 = 1 << mm2;
    tpn = 2.0 * pi / double(nn);
    mpuu2[2*k+2] = double(nn1*nn2);
    
    for (j=0;j<nn2;j++) {
      for(i=0;i<nn1;i++) {
	iu = ku + i + j * nn1;
	t1 = tpn * double(i * j);
	sin_temp = sin(t1);
	cos_temp = cos(t1);
	mpuu2[(iu-1)*2] = cos_temp;
	mpuu2[(iu-1)*2+1] = sin_temp;
      }
    }
    ku += nn;
  }

  return;
}

