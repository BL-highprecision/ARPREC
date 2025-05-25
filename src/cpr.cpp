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

int mp_real::mpcpr(const mp_real& a, const mp_real& b)
{
  /**
   * This routine compares the MP numbers A and B and returns
   * -1, 0, or 1, depending on whether A < B, A = B, or A > B
   * respectively. Debug begins with debug_level = 9.
   */
  int ia, ib;
	int i; 

  if(error_no != 0) {
    if(error_no == 99)
      mpabrt();
    return 0;
  }
  if(debug_level >= 9) {
    cerr <<"Comparing mp_real..." << endl;
  }
  ia = a[1] == 0.0 ? 0 : int(sign(1.0, a[1]));
  ib = b[1] == 0.0 ? 0 : int(sign(1.0, b[1]));
  
  if(ia != ib) 
    return ia > ib ? 1 : -1;
  
  // The signs are the same. Compare exponents.
  double ma = a[2], mb = b[2]; // grab exponents
  if(ma != mb) 
    return int(ia * sign(1.0, ma - mb));

  // The signs and exponents are the same.  Compare mantissas
  int na = std::min(int(std::abs(a[1])), prec_words);
  int nb = std::min(int(std::abs(b[1])), prec_words);

  // compare all the mantissa words that are non-zero in each variable. 
  for(i=FST_M; i< std::min(na, nb) + FST_M; i++) {
    if(a[i] != b[i])
      return ia * int(sign(1.0, a[i] - b[i]));
  }
  
  // The mantissas are the same to the common length.
  //if(na != nb)
  //  return ia * sign(1, na - nb);
	// check that trailing mantissa in longer word is zero to working precision
	if(na > nb){
		for(i=nb + FST_M; i < na + FST_M; i++){
			// if a non-zero mantissa word is found to current working precision
			// then a has a larger absolute value. return sign of a
			if( a[i] != 0.0 ) return ia;
		}
	}
	else if( nb > na ){
		for(i=na + FST_M; i < nb + FST_M; i++){
			if( b[i] != 0.0 ) return ib; 
		}
	}

  //else they are the same!
  return 0;
}

