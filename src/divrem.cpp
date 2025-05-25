#include <arprec/mp_int.h>
#include "small_inline.h"

void divrem(const mp_int &dividend, const mp_int &divisor,
		        mp_int& quotient, mp_int& remainder) {
  int nws = quotient.prec_words;
  quotient.prec_words = int(dividend[2] - divisor[2]) + 3;
  int prec_words = mp::prec_words;

  //Check to see if we have enough precision availible.
  if(quotient.prec_words > quotient.n_mantissa+6) {
    mp_int::mp_int_prec_error("divrem");
  }
  mp_real t(0.0, static_cast<int>(dividend[0]));
  mp_real::mpdivx(dividend, divisor, quotient, prec_words);
  mp_real::mpinfr(quotient, quotient, t, prec_words, 0);
  quotient.prec_words = nws;

  mp_real::mpmulx(quotient, divisor, t, prec_words);
  mp_real::mpsub(dividend, t, remainder, prec_words);

  //check to see if quotient off by one.
  //This can happen when the division should have 
  //resulted in an integer, but ended up an epsilon off,
  //or when the division was just under and integer, 
  //such as 3.999999, and the division rounded up.
  if(dividend[1] > 0.0) {
    if(remainder[1] < 0.0) {
      if(quotient[1] >= 0.0) {
	quotient--;
	mp_real::mpadd(remainder, divisor, remainder, prec_words);
      } else {
	//divisor < 0
	quotient++;
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
      }
    } else if(divisor[1] > 0.0) {
      //quotient positive
      if(remainder >= divisor) {
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
	quotient++;
      }
    } else {
      //divisor < 0, quotient negative
      remainder[1] = -remainder[1];
      if(remainder <= divisor) {
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
	quotient--;
      }
      remainder[1] = -remainder[1];
    }
  } else {
    //dividend < 0.0
    if(remainder[1] > 0.0) {
      if(quotient[1] >= 0.0) {
	//divisor < 0
	quotient--;
	mp_real::mpadd(remainder, divisor, remainder, prec_words);
      } else {
	//divisor > 0
	quotient++;
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
      }
    } else if(divisor[1] < 0.0) {
      //quotient positive
      if(remainder <= divisor) {
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
	quotient++;
      }
    } else {
      //divisor > 0, quotient negative
      remainder[1] = -remainder[1];
      if(remainder >= divisor) {
	mp_real::mpsub(remainder, divisor, remainder, prec_words);
	quotient--;
      }
      remainder[1] = -remainder[1];
    }
  }

  return;
}

