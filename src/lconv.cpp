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

void mp_real::mplconv(int iq, int n, int nsq, double *a, double *b, double *c)
{
  /**
   * This computes the linear convolution of the N-long inputs A and B.
   * |IQ| is the number of arguments (i.e., if IQ = 1, then B is ignored).
   * If IQ is negative (and N < 64) then only the second half of the result
   * vector is required (i.e. this is a call by itself -- see below).
   * NSQ = int(sqrt(3*N)) is an input required for the dimention of
   * DC1 and DC2 (see below).
   *
   * This routine employs an advanced FFT-based scheme, except for
   * small n.  This rouine is not intende to be called directly by the user.

   * Two machine-dependent parameters are ste in this routine:
   *    ERM = Maximum tolerated FFT roundoff erorr.  On IEEE systems ERM = 
   *    0.438.  It is not necessary to specify ERM for modest levels of
   *    precision -- see comments below.
   *    MBT = Number of mantissa bits in double precision data.  MBT = 53
   *    since this library is intended for IEEE systems only (without
   *    extended registers being used).
   */
  const double cl2 = 1.4426950408889633;
  int ncr1, k, n1, n2, m, m2, m21, m1, nm, i, j, n4;
  double t1 ;
  double td1;
  double temp;

  // Handle the case where N is less than ncr1 = pow(2, mpmcrx-1).
  // if IQ < 0, only the second half of the result vector is returned,
  // since the first half won't be used.

  ncr1 = 1 << (mpmcrx-1);
  if(n < ncr1) {
    if (iq == 1) {
      for(k=0;k<2*n;k++) {
	td1 = 0.0;
	n1 = std::max (k - n +1, 0);
	n2 = std::min(k+1, n);
	
	for(j=n1;j<n2;j++)
	  td1 += a[j] * a[k-j];

	c[k] = td1;
      }
    } else if(iq == 2) {
      for(k=0; k < 2*n; k++) {
	td1 = 0.0;
	n1 = std::max (k - n +1, 0);
	n2 = std::min(k+1, n);
	
	for(j=n1;j<n2;j++)
	  td1 += a[j] * b[k-j];

	c[k] = td1;
      }
    } else if (iq == -1) {
      for(k=0;k<n-1;k++)
	c[k] = 0.0;

      for(k=n-1; k < 2*n; k++) {
	td1 = 0.0;
	n1 = std::max (k - n +1, 0);
	//n2 == n always

	for(j=n1;j<n;j++)
	  td1 += a[j] * a[k-j];

	c[k] = td1;
      }
    } else if (iq == -2) {
      for(k=0;k<n-1;k++)
	c[k] = 0.0;

      for(k=n-1; k < 2*n; k++) {
	td1 = 0.0;
	n1 = std::max (k - n +1, 0);
	//n2 == n always

	for(j=n1;j<n;j++)
	  td1 += a[j] * b[k-j];

	c[k] = td1;
      }
    }
    return;
  }
  
  double *dc1 = new double[6*n + 2*nsq*mpnsp1+6];
  double *dc2 = new double[6*n + 2*nsq*mpnsp1+6]; //for complex numbers
  double *d1 = new double[3*n+2];
  double *d2 = new double[3*n+2];
  double *d3 = new double[3*n+2];


  // Determine M1 and N1.  Note that by this reckoning, N1 <= 1.5 N.
  // This is the reason for the 3*n/2 dimensions above.

  t1 = 0.75 * n;
  m1 = int(cl2 * log(t1) + 1 - mprxx);
  n1 = 1 << m1;
  m2 = m1 + 1;
  n2 = 2 * n1;
  n4 = 2 * n2;
  nm = std::min(2 * n, n2);
  if(std::abs(iq) == 1) {
    for(i=0; i<n; i++)
      d1[i] = a[i];

    for(; i<n1; i++)
      d1[i] = 0.0;

    // Perform a forward real-to complex FFT on the vector in a.
    
    mpfftrc(1, m2, n2, n, d1, dc1);

    // square the resulting complex vector;
    for(i=0;i< 2*(n1+1);i+=2) {
      temp = dc1[i] * dc1[i+1];
      dc1[i]  =  sqr( dc1[i] ) - sqr ( dc1[i+1] );
      dc1[i+1] = 2.0 * temp;
    }
      
  } else {
    //std::abs(iq) == 2.
    for(i=0;i<n;i++) {
      d1[i] = a[i];
      d2[i] = b[i];
    }
    
    for(i=n;i<n1;i++)
      d1[i] = d2[i] = 0.0;


    // Perform forward real-to-complex FFTs on the vectors in a and in b.
    mpfftrc(1, m2, n2, n, d1, dc1);
    mpfftrc(1, m2, n2, n, d2, dc2);
    

    // Multiply the resulting complex vectors.
      for(i=0;i<(n1+1)*2;i+=2) {
        temp = dc1[i] * dc2[i] - dc1[i+1] * dc2[i+1];
        dc1[i+1] = dc1[i] * dc2[i+1] + dc1[i+1] * dc2[i];
        dc1[i] = temp;
      }
  }

  // Perform an inverse complex-to-real FFT on the resulting data.
  
  mpfftcr(-1, m2, n2, nsq, dc1, dc2);

  // divide by N4.
  double an = 1.0 / n4;
  int ms;
  double td2;
 
  //usually this line is commented out. see below. 
  //#define MPMULX_ROUNDOFF_ERROR_TEST

  for(i=0;i<nm;i++) {
    td1 = an * (dc2[i]);
    td2 = nint(td1); 
#ifdef MPMULX_ROUNDOFF_ERROR_TEST
    d1[i] = std::abs(td2 - td1);
#endif
    c[i] = td2;
  }
  
  /**
   * Find the largest FFT roundoff error.  Roundoff error is minimal unless
   * exceedingly high precision (i.e. over one million digits) is used.  Thus
   * this test may be disabled in normal use.  To disable this test,
   * comment out the above definition of MPMULX_ROUNDOFF_ERROR_TEST.
   * 
   * This code can be used as a rigorous system integrity test.
   * First set MBT according to the system being used (see above),
   * then set erm to be fairly small, say 0.001 or whatever is somewhat
   * larger that the largest FFT roundoff error typically encountered for 
   * a given precision level on the computer being used.  Enable this
   * test as explained in the previous paragraph.  Then if an anomalously 
   * large roundoff error is detected, a ahardware or compiler error
   * has likely occurred.
   */
#ifdef MPMULX_ROUNDOFF_ERROR_TEST
  {
    const double erm = 0.438, mbt = 53;
    int i2, i3, i4, i5;
    
    t1 = 0.0;

    for(i=0;i<nm;i++) {
      if(d1[i] > t1) {
	i1 = i;
	t1 = d1[i];
      }
    }
  
    // Check if maximum roundoff error exceeds the limit erm, which is set
    // above.  Also determine the number of fractional bits and how large
    // the error is in terms of units in the last place (ulp).
    if(t1 > erm) {
      if(MPKER[55] != 0) {
	t2 = an * d1[i1];
	i2 = cl2 * log(t1) + 1.0 + mprxx;
	i3 = cl2 * log(t2) + 1.0 + mprxx;
	i4 = mbt + i2 - i3;
	i5 = t1 * (1<<i4) + mprxx;
	cerr << "*** MPLCONV : Excessive FFT roundoff error: "
	     <<i1 <<" "<<t1<<" "<<i4<<" "<<i5 << endl;
	error_no = 55;
	if(MPKER[error_no] == 2) mpabrt();
      }
    } 
  } 
#endif
 
  // Handle case where n > n1;
  // In other words, fix a few words if the convolution if
  // It was decided above that it was better to do an
  // FFT of size smaller than n. This can happen
  // when n is just a few more than a power of two.
  // The FFTs are done on power of two sized arrays only.

  if(n > n1) {
    m = n - n1;
    m2 = 2 * m;
    m21 = 2 * m - 1;
    ms = int(ANINT(sqrt(3.0 * m21) + mprxx));
    k = n1 - m + 1;
    if(std::abs(iq) == 1) {
      for(i=0;i<m21;i++) 
	d1[i] = a[k+i];
      
      mplconv(-1, m21, ms, d1, d2, d3);
    } else {
      for(i=0;i<m21;i++) {
	d1[i] = a[k+i];
	d2[i] = b[k+i];
      }
      mplconv(-2, m21, ms, d1, d2, d3);
    }
    int ii;
    for(i=0;i<m2;i++) {
      ii = i + m2 - 2;
      c[i] -= d3[ii];
      c[i+n2] = d3[ii];
    }
  }
  delete [] dc1;
  delete [] dc2;
  delete [] d1;
  delete [] d2;
  delete [] d3;

  return;
} // end mp_real::mplconv

