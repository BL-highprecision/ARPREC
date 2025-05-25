/*
 * include/arprec/mp_int.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef __MPINT_H
#define __MPINT_H

#include <arprec/mp_real.h>
#include <arprec/mp_int_temp.h>

/**
 * This class representes MP integer numbers.
 */

struct ARPREC_API mp_int: public mp_real {

  //destructor is unneccessary. The base class will do all the work.

  //Constructors 

  mp_int():mp_real() {}

  //This constructor is different than mp_real!!!!
  mp_int(int s):mp_real() {
    if(s > 0) {
      mpr[1] = 1.0;
      mpr[2] = 0.0;
      mpr[3] = double(s);
    } else if(s == 0) {
      zero(*this);
    } else {
      mpr[1] = -1.0;
      mpr[2] = 0.0;
      mpr[3] = double(-s);
    }
  }
   mp_int(double d):mp_real(aint(d)) {};
   mp_int(const mp_real& a):mp_real(0.0, static_cast<int>(a[0]))
     {
       if(std::abs(a[1]) > a[2]+1.0) {
         int nws = prec_words;
         prec_words = std::min(int(a[2])+1, prec_words);
         if (prec_words < 1)
          prec_words = 1;
         mp_real t;
         mpinfr(a, *this, t, prec_words);
         prec_words = nws;
      } else {
	mpeq(a, *this, prec_words);
      }
    }
  mp_int(const mp_int& a):mp_real(a) {}
  mp_int(const mp_int_temp& other):mp_real(other.to_mp_real_temp()) {}
  mp_int(double *dp):mp_real(dp) {
    if(mpr[1] != 0.0) {
      if(std::abs(mpr[1]) > mpr[2]+1.0) {
	//truncate
	if(mpr[1] > 0.0) { 
	  mpr[1] = mpr[2] + 1.0;
	} else {
	  mpr[1] = -mpr[2] -1.0;
	}
      }
    }
  }
  mp_int(const double *dp):mp_real(dp) {
    if(mpr[1] != 0.0) {
      if(std::abs(mpr[1]) > mpr[2]+1.0) {
	//truncate
	if(mpr[1] > 0.0) { 
	  mpr[1] = mpr[2] + 1.0;
	} else {
	  mpr[1] = -mpr[2] -1.0;
	}
      }
    }
  }

  //Assignment operators.

  mp_int &operator=(const mp_int &a) {
    mpeq(a, *this, prec_words);
    mp_int::ovcheck(*this);
    return *this;
  }
  mp_int &operator=(const mp_real_temp& a) {
    //use superclass version.
    mp_real::operator=(a);
    if(mpr[1] > mpr[2]+1.0)
      mpr[1] = mpr[2]+1.0;
    mp_int::ovcheck(*this);
    return *this;
  }
  mp_int &operator=(const mp_int_temp& a) {
    //use superclass version.
    mp_real::operator=(a.to_mp_real_temp());
    mp_int::ovcheck(*this);
    return *this;
  }
  mp_int &operator=(int a) {
    if(a == 0) {
      mp_real::zero(*this);
    } else if(a > 0) {
      mpr[1] = 1.0; mpr[2] = 0.0; mpr[3] = static_cast<double>(a);
    } else {
      mpr[1] = -1.0; mpr[2] = 0.0; mpr[3] = -static_cast<double>(a);
    }
    mp_int::ovcheck(*this);
    return *this;
  }

  mp_int_temp toTempAndDestroy() { 
    mp_int_temp ret(this->mpr); mpr = 0; return ret;
  }
  
  
  /**
   * integer problems require exactness, without rounding.
   * if not enough precision is availible, this will be called.
   */
  static void mp_int_prec_error(const char *routine=static_cast<const char *>("")) {
    std::cerr << "\n*** MPINT ERROR: not enough precision for "<< routine<<".";
    mpabrt();
  }

  /*
   * Check integer overflow.
   */
  static void ovcheck(const mp_int &ia) {
    if ((ia[0]-5) <= ia[2]) {
      if (MPKER[74] != 0) {
        std::cout << "\n*** MP integer overflow.";
	error_no = 74;
	if (MPKER[error_no] == 2) mpabrt();
      }
    }
    return;
  }

  /* divides, produces integer quotient
   * and also returns the remainder */
  
  mp_int& operator++();
  mp_int& operator++(int);
  mp_int& operator--();
  mp_int& operator--(int);


  mp_int& operator+=(const mp_int& other) {
    mp_real::operator+=(other);
    return *this;
  }
  mp_int& operator+=(int other) {
    mp_real::operator+=(mp_int(other));
    return *this;
  }
  mp_int& operator-=(int other) {
    mp_real::operator-=(mp_int(other));
    return *this;
  }
  mp_int& operator*=(int other) {
    mp_real::operator*=(double(other));
    return *this;
  }

  mp_int& operator/=(const mp_int& other);
  mp_int& operator/=(int other);

  mp_int& operator%=(const mp_int& other);
  mp_int& operator%=(int other);

};

ARPREC_API void divrem(const mp_int &dividend, const mp_int &divisor,
  mp_int &quotient, mp_int &remainder);
ARPREC_API mp_int_temp pow(const mp_int &a, int n);
ARPREC_API mp_int_temp pow(const mp_int &a, const mp_int &b);

ARPREC_API mp_int_temp operator+(const mp_int &a, const mp_int& b);
ARPREC_API mp_int_temp operator+(int a, const mp_int& b);
ARPREC_API mp_int_temp operator+(const mp_int& a, int b);
ARPREC_API mp_real_temp operator+(const mp_int& a, double b);
ARPREC_API mp_real_temp operator+(double a, const mp_int& b);

ARPREC_API mp_int_temp operator-(const mp_int &a, const mp_int& b);
ARPREC_API mp_int_temp operator-(int a, const mp_int& b);
ARPREC_API mp_int_temp operator-(const mp_int& a, int b);
ARPREC_API mp_real_temp operator-(const mp_int& a, double b);
ARPREC_API mp_real_temp operator-(double a, const mp_int& b);

ARPREC_API mp_int_temp operator*(const mp_int &a, const mp_int& b);
ARPREC_API mp_int_temp operator*(int a, const mp_int& b);
ARPREC_API mp_int_temp operator*(const mp_int& a, int b);
ARPREC_API mp_real_temp operator*(const mp_int& a, double b);
ARPREC_API mp_real_temp operator*(double b, const mp_int& a);

ARPREC_API mp_int_temp operator/(const mp_int &a, const mp_int& b);
ARPREC_API int operator/(int a, const mp_int& b);
ARPREC_API mp_int_temp operator/(const mp_int& a, int b);
ARPREC_API mp_real_temp operator/(const mp_int& a, double b);
ARPREC_API mp_real_temp operator/(double a, const mp_int& b);

ARPREC_API mp_int_temp operator%(const mp_int &a, const mp_int& b);
ARPREC_API int operator%(int a, const mp_int& b);
ARPREC_API int operator%(const mp_int& a, int b);

//truncate to integer. 
//This function silently returns a bad value
//if a cannot be stored in an integer.
ARPREC_API int integer(const mp_int& a);

ARPREC_API bool operator==(const mp_int& a, const mp_int& b);
ARPREC_API bool operator!=(const mp_int& a, const mp_int& b);
ARPREC_API bool operator>(const mp_int& a, const mp_int& b);
ARPREC_API bool operator<(const mp_int& a, const mp_int& b);
ARPREC_API bool operator>=(const mp_int& a, const mp_int& b);
ARPREC_API bool operator<=(const mp_int& a, const mp_int& b);

ARPREC_API bool operator==(const mp_int& a, const mp_real& b);
ARPREC_API bool operator!=(const mp_int& a, const mp_real& b);
ARPREC_API bool operator>(const mp_int& a, const mp_real& b);
ARPREC_API bool operator<(const mp_int& a, const mp_real& b);
ARPREC_API bool operator>=(const mp_int& a, const mp_real& b);
ARPREC_API bool operator<=(const mp_int& a, const mp_real& b);

ARPREC_API bool operator==(const mp_real& a, const mp_int& b);
ARPREC_API bool operator!=(const mp_real& a, const mp_int& b);
ARPREC_API bool operator>(const mp_real& a, const mp_int& b);
ARPREC_API bool operator<(const mp_real& a, const mp_int& b);
ARPREC_API bool operator>=(const mp_real& a, const mp_int& b);
ARPREC_API bool operator<=(const mp_real& a, const mp_int& b);
  
#if (ARPREC_INLINE)
#include <arprec/mp_int_inline.h>
#else
bool operator==(const mp_int& a, const mp_real& b);
bool operator!=(const mp_int& a, const mp_real& b);
bool operator>(const mp_int& a, const mp_real& b);
bool operator<(const mp_int& a, const mp_real& b);
bool operator>=(const mp_int& a, const mp_real& b);
bool operator<=(const mp_int& a, const mp_real& b);
bool operator==(const mp_real& a, const mp_int& b);
bool operator!=(const mp_real& a, const mp_int& b);
bool operator>(const mp_real& a, const mp_int& b);
bool operator<(const mp_real& a, const mp_int& b);
bool operator>=(const mp_real& a, const mp_int& b);
bool operator<=(const mp_real& a, const mp_int& b);
bool operator==(const mp_int& a, const mp_int& b);
bool operator!=(const mp_int& a, const mp_int& b);
bool operator>(const mp_int& a, const mp_int& b);
bool operator<(const mp_int& a, const mp_int& b);
bool operator>=(const mp_int& a, const mp_int& b);
bool operator<=(const mp_int& a, const mp_int& b);
void divrem(const mp_int &dividend, const mp_int &divisor, mp_int& quotient, mp_int& remainder);
mp_int_temp operator+(const mp_int& a, const mp_int& b);
mp_int_temp operator+(const mp_int& a, int b);
mp_int_temp operator+(int b, const mp_int& a);
mp_real_temp operator+(const mp_int& a, double b);
mp_real_temp operator+(double b, const mp_int& a);
mp_int_temp operator-(const mp_int& a, const mp_int& b);
mp_int_temp operator-(const mp_int& a, int b);
mp_int_temp operator-(int b, const mp_int& a);
mp_real_temp operator-(const mp_int& a, double b);
mp_real_temp operator-(double b, const mp_int& a);
mp_int_temp operator*(const mp_int& a, const mp_int& b);
mp_int_temp operator*(const mp_int& a, int b);
mp_int_temp operator*(int b, const mp_int& a);
mp_real_temp operator*(const mp_int& a, double b);
mp_real_temp operator*(double b, const mp_int& a);
mp_int_temp operator/(const mp_int& a, const mp_int& b);
mp_int_temp operator/(const mp_int& a, int b);
int operator/(int a, const mp_int& b);
mp_real_temp operator/(const mp_int& a, double b);
mp_real_temp operator/(double b, const mp_int& a);
mp_int_temp operator%(const mp_int& a, const mp_int &b);
int operator%(int a, const mp_int &b);
int operator%(const mp_int& a, int b);
mp_int_temp pow(const mp_int &a, int n);
mp_int_temp pow(const mp_int &a, const mp_int& b);
int integer(const mp_int& a);

#endif

#endif
