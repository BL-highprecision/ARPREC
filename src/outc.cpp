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
#define _CRT_SECURE_NO_DEPRECATE
#include <cassert>
#include <arprec/mp_real.h>
#include "small_inline.h"
#include <cstdio>

using std::cerr;
using std::endl;

/* This routine is deprecated.  
 * Use more C++ friendly mp_real::to_string or mp_real::write.  */

void mp_real::mpoutc(const mp_real& a, char* b, int& n, int prec_words)
{
  /**
   * Converts the MP number A into character form in the char array B.
   * N (an output parameter) is the length of the output.  In other words,
   * A is contained in B(1), ..., B(N).  The format is analogous to the
   * Fortran exponential format (E format), except that the exponent is
   * placed first.
   * Debug output starts with debug_level = 7.
   *
   * Max byte space for B: 15.05 * MPNW + 30 cells.
   *                              ( log10(2^50) * prec_words + 30 )
   * This routine is called by MPOUT, but it may be directly called by the
   * user if desired for custom output.
   *
   */
  int i, j, nn, no, nws;
  const double al2 = 0.301029995663981195; // log10(2)
  double aa, t1;
  char* ca = new char[17];
  int n6 = prec_words + 8;
  mp_real f(10.0, 9), sk0(0.0, n6), sk1(0.0, n6);
  int BreakLoop, Loop;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    b[0] = '\0';
    n = 0;
    return;
  }
	
  if (debug_level >= 7) print_mpreal("MPOUTC I ", a);
	
  int ia = int(sign(1.0, a[1]));
  int na = std::min (int(std::abs(a[1])), prec_words);
  nws = prec_words;
  prec_words = prec_words + 1;
  
  //  Determine exact power of ten for exponent.
  int nx;
  if (na != 0) {
    aa = a[3];
    if (na >= 2) aa = aa + mprdx * a[4];
    if (na >= 3) aa = aa + mprx2 * a[5];
    if (na >= 4) aa = aa + mprdx * mprx2 * a[6];
    t1 = al2 * mpnbt * a[2] + log10 (aa);
    if (t1 >= 0.0) 
      nx = int(t1); // *cast*
    else 
      nx = int(t1 - 1.0); // *cast*
#if 0
    printf("mpoutc[1]: aa = %22.18e, t1 = %22.18e, nx = %d\n", aa, t1, nx);
#endif
    if(nx >=0) {
      mpnpwr(f, nx, sk0, prec_words);
      mpdiv(a, sk0, sk1, prec_words);
    } else {
      mpnpwr(f, -nx, sk0, prec_words);
      mpmul(a, sk0, sk1, prec_words);
    }
    // If we didn't quite get it exactly right, multiply or divide by 10 to fix.
#if 0
    printf("mpoutc[2]: sk0, sk1\n");
    print_mpreal("  mpoutc", sk0);
    print_mpreal("  mpoutc", sk1);
#endif
    Loop = 1;
    while (Loop) {
      if (sk1[2] < 0) { // exponent < 0
        --nx;
        mpmuld(sk1, 10.0, 0, sk0, prec_words);
        mpeq(sk0, sk1, prec_words);
      } else if (sk1[3] >= 10.) {
        ++nx;
        mpdivd (sk1, 10.0, 0, sk0, prec_words);
        mpeq (sk0, sk1, prec_words);
      } else {
	Loop = 0;
      }
    }
    sk1[1] = std::abs(sk1[1]);
  } else {
    nx = 0;
  }
#if 0
  printf("mpoutc[3]: sk1\n");
  print_mpreal("  mpoutc ", sk1);
#endif
  // Now sk1 < 10
  // Place exponent first instead of at the very end as in Fortran.
  b[0] = '1';
  b[1] = '0';
  b[2] = ' ';
  b[3] = '^';
  sprintf(ca, "%10d", nx);
  int len = static_cast<int>(strlen(ca));
  int blank = 14-len;
  for(i = 4; i < blank; i++) b[i]=' ';
  for(i = 0; i < len; i++) b[blank+i] = ca[i]; 
  b[14] = ' ';
  b[15] = 'x';
  b[16] = ' ';
	
  //  Insert sign and first digit.
  if (ia == -1) b[17] = '-';
  else b[17] = ' ';
  if (na != 0) 
    nn = int(sk1[3]); // in [1, 10)
  else nn = 0;
  sprintf(ca, "%1d", nn);   
  b[18] = ca[0];
  b[19] = '.';
  int ix = 20;
  if (na == 0) {     
    n = ix;
    prec_words = nws;
    if (debug_level >= 7) {
      no = std::min (n, 6 * debug_words + 20);
      cerr << "MPOUTC O "; 
      for(i = 0; i < no; i++) cerr << b[i];
      cerr << endl;
    }
    delete [] ca;
    b[n]='\0';
    return;
  }
	
  f[3] = double(nn); // *cast*
  mpsub (sk1, f, sk0, prec_words); // remainder
#if 0 
  printf("mpoutc[4]: f, sk0\n");
  print_mpreal("  mpoutc[4] ", f);
  print_mpreal("  mpoutc[4] ", sk0);
#endif
  if (sk0[1] == 0) {     
    n = ix;
    if (debug_level >= 7) {
      no = std::min (n, 6 * debug_words + 20);
      cerr << "MPOUTC O "; 
      for(i=0; i<no; i++)
        cerr << b[i];
      cerr << endl;
    }
    delete [] ca;
    b[n]='\0';
    return;
  }

  mpmuld (sk0, 1e6, 0, sk1, prec_words);
  int nl = int(std::max (prec_words * log10 (mpbdx) / 6.0 - 1.0, 1.0));
  int mpnw_change = int(mpnbt * log(double(2.0))/log(double(10.0))/6.0)+2;
#if 0
  printf("mpoutc[5]: mpnw_change %d, sk1 = \n", mpnw_change);
  print_mpreal("  mpoutc[5] ", sk1);
#endif
  //  Insert the digits of the remaining words. 6 decimal digits at a time.
  BreakLoop = 0;
  for (j = 1; j <= nl; ++j) {
    if (sk1[2] == 0.) {
      assert(sk1[3] <= 2e9);
      nn = int(sk1[3]);
      f[1] = 1;
      f[3] = double(nn);
    } else { // exponent < 0
      f[1] = 0;
      nn = 0;
    }
    
    sprintf(ca, "%6d", nn);
#if 0
    printf("mpoutc[6]: %6d\n", nn);
#endif
    for (i = 0; i < 6; ++i) b[i+ix] = (ca[i] != ' ' ? ca[i] : '0');
    
    ix += 6;
    mpsub(sk1, f, sk0, prec_words);
    mpmuld(sk0, 1e6, 0, sk1, prec_words);
    if (sk1[1] == 0) {
      BreakLoop = 1;
      break;
    }
    if(!((j) % (mpnw_change))) {
      prec_words--;
    }
  }

  //  Check if trailing zeroes should be trimmed.
  if (!BreakLoop) j = nl + 1;
  
  int l = --ix;
  if (b[l] == '0' ||
      (j > nl && b[l-1] == '0' && b[l-2] == '0' && b[l-3] == '0')) {
    b[l] = '\0';
    BreakLoop = 0;
    for (i = l - 1; i >= 20; --i) {
      if (b[i] != '0') {
        ix = i;
	BreakLoop = 1;
        break;
      }
      b[i] = '\0';
    }

    if (!BreakLoop) ix = 20;
    
  } else if (j > nl && b[l-1] == '9' && b[l-2] == '9' && b[l-3] == '9') {
    // Check if trailing nines should be rounded up.
    b[l] = '\0';
    int roundUp=0;

    BreakLoop = 0;
    for (i = l - 1; i >= 20; --i) {
      if (b[i] != '9') {
        BreakLoop = 1;
        break;
      }
      roundUp = 1;
      b[i] = '\0';
    }
    
    // We have rounded away all digits to the right of the decimal point,
    // and the digit to the left of the digit is a 9.  Set the digit to 1
    // and increase the exponent by one.
    if (!BreakLoop) {
      ix = 20;
      if (b[18] == '9') {
        b[18] = '1';
        sprintf(ca, "%10d", nx+1);
        for (i = 0; i < 10; ++i) b[i+4] = ca[i];   
      } else {
        ca[0] = b[18]; ca[1]='\0';
        nn = atoi(ca);
        sprintf(ca, "%1d", nn+1);
        b[18] = ca[0];
      }
    } else if(roundUp){
      ca[0] = b[i]; ca[1]='\0';
      nn = atoi(ca);
      sprintf(ca, "%1d", nn+1);
      b[i] = ca[0];
      ix = i;
    }
  } 
  
  n = ix;
  if (debug_level >= 7) {
    no = std::min (n, 6 * debug_words + 20);
    cerr << "MPOUTC O "; 
    for(i = 0; i <= no; ++i) cerr << b[i];
    cerr << endl;
  }
  b[++n]='\0';
  prec_words = nws;
  delete [] ca;
  return;
}

