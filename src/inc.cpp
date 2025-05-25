#include <arprec/mp_int.h>
#include "small_inline.h"

mp_int& mp_int::operator++()
{
  int prec_words = mp::prec_words;

  mp_int *a = new mp_int(*this);
  if(mpr[1] == 0.0) {
    mpr[1] = 1.0;
    mpr[2] = 0.0;
    mpr[FST_M] = 1.0;
    return *a;
  }

  if(mpr[2]+1.0 < std::abs(mpr[1])) {
    mp_int_prec_error("operator++()");
  }
  if(mpr[2]+1.0 > std::abs(mpr[1])) {
    mpadd(*this, mp_real(1.0), *this, prec_words);
    return *a;
  } else {
    double b;
    if(mpr[1] >= 0) {
      b = (mpr[FST_M + int(std::abs(mpr[1])) - 1] += 1.0);
      if(b >= mpbdx) {
	//Take the pointer of the function mpnorm in the hope that
	//it will not expand inline here
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *a;
    } else {
      b = (mpr[FST_M + int(std::abs(mpr[1])) - 1] -= 1.0);
      if(b <= 0.0) {
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *a;
    }
  }
}

mp_int& mp_int::operator++(int)
{
  int prec_words = mp::prec_words;

  if(mpr[1] == 0.0) {
    mpr[1] = 1.0;
    mpr[2] = 0.0;
    mpr[FST_M] = 1.0;
    return *this;
  }

  if(mpr[2]+1.0 < std::abs(mpr[1])) {
    mp_int_prec_error("operator++(int)");
  }
  if(mpr[2]+1.0 > std::abs(mpr[1])) {
    mpadd(*this, mp_real(1.0), *this, prec_words);
    return *this;
  } else {
    double b;
    if(mpr[1] >= 0) {
      b = (mpr[FST_M + int(std::abs(mpr[1])) - 1] += 1.0);
      if(b >= mpbdx) {
	//Take the pointer of the function mpnorm in the hope that
	//it will not expand inline here
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *this;
    } else {
      b = (mpr[FST_M + int(std::abs(mpr[1])) - 1] -= 1.0);
      if(b <= 0.0) {
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *this;
    }
  }
}

mp_int& mp_int::operator--()
{
  mp_int *a = new mp_int(*this);
  int prec_words = mp::prec_words;

  if(mpr[1] == 0.0) {
    mpr[1] = -1.0;
    mpr[2] = 0.0;
    mpr[FST_M] = 1.0;
    return *a;
  }

  if(mpr[2]+1.0 < std::abs(mpr[1])) {
    mp_int_prec_error("operator--()");
  }
  if(mpr[2]+1.0 > std::abs(mpr[1])) {
    mpadd(*this, mp_real(-1.0), *this, prec_words);
    return *a;
  } else {
    double b;
    if(mpr[1] < 0) {
      b = (mpr[FST_M + int(std::abs(mpr[1]))] += 1.0);
      if(b >= mpbdx) {
	//Take the pointer of the function mpnorm in the hope that
	//it will not expand inline here
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *a;
    } else {
      b = (mpr[FST_M + int(std::abs(mpr[1]))] -= 1.0);
      if(b <= 0.0) {
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *a;
    }
  }
}

mp_int& mp_int::operator--(int)
{
  int prec_words = mp::prec_words;

  if(mpr[1] == 0.0) {
    mpr[1] = -1.0;
    mpr[2] = 0.0;
    mpr[FST_M] = 1.0;
    return *this;
  }

  if(mpr[2]+1.0 < std::abs(mpr[1])) {
    mp_int_prec_error("operator--(int)");
  }
  if(mpr[2]+1.0 > std::abs(mpr[1])) {
    mpadd(*this, mp_real(-1.0), *this, prec_words);
    return *this;
  } else {
    double b;
    if(mpr[1] < 0) {
      b = (mpr[FST_M + int(std::abs(mpr[1]))] += 1.0);
      if(b >= mpbdx) {
	//Take the pointer of the function mpnorm in the hope that
	//it will not expand inline here
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *this;
    } else {
      b = (mpr[FST_M + int(std::abs(mpr[1]))] -= 1.0);
      if(b <= 0.0) {
	(&mpnorm)(mpr, *this, prec_words);
      }
      return *this;
    }
  }
}

