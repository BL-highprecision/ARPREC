#include <arprec/mp_int.h>
#include "small_inline.h"

mp_int_temp operator/(const mp_int& a, const mp_int& b) {
  mp_int c, d;
  int prec_words = mp::prec_words;

  if(a[1] == 0.0) {
    mp_real::zero(c);
    return c.toTempAndDestroy();
  }

  int nws = c.prec_words;
  c.prec_words = int(a[2] - b[2]) + 4;
  d.prec_words = int(a[2] - b[2]) + 4;
  if(c.prec_words > c.n_mantissa+1) {
    mp_int::mp_int_prec_error("operator/(mp_int, mp_int)");
  }
  mp_real t;
  enum mp::rounding_mode roun = mp::round_dir;
  mp_real::round_dir = mp::round_to_zero; //set to round toward zero. works most of the time.
  mp_real::mpdivx(a, b, d, prec_words);
  mp_real::mpinfr(d, c, t, prec_words);
  mp_real::round_dir = roun; //restore round_dir.

  //Check to see if one off.
  //A check is needed in exact and almost exact cases.
  t[1] = std::abs(t[1]);
  if(t[1] == 0.0 || t > 0.999) {
    c.prec_words = nws;
    mp_int t2;
    mp_int::mpmulx(c, b, t2, prec_words);
    mp_int::mpsub(a, t2, t2, prec_words);
    if(a[1] > 0.0) {
      if(b[1] > 0.0) {
	if(t2 < mp_int(0)) {
	  c--;
	} else if(t2 >= b) {
	  c++;
	} 
      } else {
	if(t2 < mp_int(0)) {
	  c++;
	} else if(t2 >= -b) {
	  c--;
	}
      }
    } else {
      if(b[1] > 0.0) {
	if(t2 > mp_int(0)) {
	  c++;
	} else if(t2 <= -b) {
	  c--;
	}
      } else {
	if(t2 > mp_int(0)) {
	  c--;
	} else if(t2 <= b) {
	  c++;
	}
      }
    }
  }

  c.prec_words = nws;
  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

mp_int_temp operator/(const mp_int& a, int b) {
  mp_int c(a[0]);
  int nws = c.prec_words;
  int prec_words = mp::prec_words;

  c.prec_words = int(a[2]) + 3;
  if(c.prec_words > c.n_mantissa+6) {
    mp_int::mp_int_prec_error("operator/(mp_int, int)");
  }

  mp_real t;
  mp_real::mpdivd(a, b, 0, c, prec_words);
  mp_real::mpinfr(c, c, t, prec_words, 0);
  c.prec_words = nws;

  mp_int::ovcheck(c);
  return c.toTempAndDestroy();
}

int operator/(int a, const mp_int& b) {
  if(b[1] == 0.0) {
    std::cerr << "\n*** MPINT, operator/(int, mp_int) : division by zero";
    mp_int::mpabrt();
  }
  if(b[2] > 0.0 || b[FST_M] > std::abs(double(a))) {
    return 0 ;
  }
  return a / int(b[FST_M] * b[1]);
}

