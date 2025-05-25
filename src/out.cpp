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

void mp_real::mpout(const mp_real& a, int la, char* cs, int& l, int prec_words)
{
  /**
   * This routine writes the exponent plus LA mantissa digits of the MP number
   * A to the character array CS. CS must
   * be dimensioned at least LA + 25.  The digits of A may span more than one
   * line (newlines are included in the array).  
   * A comma is placed at the end of the last line to denote the end of
   * the MP number.  Here is an example of the output:
   *
   * 10 ^  -4 x  3.14159265358979323846264338327950288419716939937510,
   */

  int ll, nws;
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    return;
  }
  
  nws = prec_words;
  ll = int(la / log10 (mpbdx) + 3.0);
  prec_words = std::min (prec_words, ll);
#if 0
  printf("mpout[1] prec_words %d\n", prec_words); 
#endif
  mpoutc(a, cs, l, prec_words);
#if 0
  printf("mpout[2] after mpoutc l %d\n", l);
  {
    int i;
    for (i = 0; i < l; ++i) printf("%c", cs[i]);
    printf("\n");
  }
#endif
  prec_words = nws;
  l = std::min (l, la + 20) + 1;
  if(isdigit(cs[l-1]) && cs[l-1] >= '5') {
    int i;
    for(i=l-2;i>0 && cs[i] == '9'; i--);
    if(!i) {
      //This shouldn't happen.
      assert(0);
    }
    if(cs[i] == '.') {
      if(cs[i-1] == '9') {
	//This is a bug.  This should push the exponent up by one.
	//instead, it just leaves a bunch of nines.
	long exponent;
	int n=i, j= l-1;
	for(i=0;i<j && cs[i] != '^';i++);
	if(i == j){
	  cerr << "*** MPOUT : error when rounding 9's" << endl;
	  cs[l-1] = ',';
	  cs[l] = '\0';
	  return;
	}
	for(;j>i && cs[j] != 'x';j--);
	if(i == j) {
	  cerr << "*** MPOUT : error when rounding 9's" << endl;
	  cs[l-1] = ',';
	  cs[l] = '\0';
	  return;
	}
	cs[j] = '\0';
	exponent = atol(cs+i+1);
	cs[j] = 'x';
	exponent++;
	sprintf(cs + (j-10), "%10ld", exponent);
	
	cs[n-1] = '1';
	cs[n+1] = ',';
	cs[n+2] = '\0';
	l = n+2;
      } else {
	cs[i-1]++;
	cs[i+1] = ',';
	cs[i+2] = '\0';
	l = i+2;
      }
    } else {
      cs[i]++;
      cs[i+1] = ',';
      cs[i+2] = '\0';
      l = i+2;
    }
  } else {
    int i;
    for(i=l-2;i>0 && cs[i] == '0'; i--);
    cs[i+1] = ',';
    cs[i+2] = '\0';
    l = i+2;
  }
  return;
}

