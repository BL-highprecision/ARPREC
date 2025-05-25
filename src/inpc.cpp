/*
 * src/mpreal4.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 * Additional routines, mostly those for which speed is less important.
 */
#include <cstdlib>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

void mp_real::mpinpc (const char a[], int n, mp_real& b)
{
  /**
   * Converts the char[] array A of length N into the MP number B.  The
   * string A must be in the format '10^s a x tb.c' where a, b and c are digit
   * strings; s and t are '-', '+' or blank; x is either 'x' or '*'.  Blanks
   * may be embedded anywhere.  The exponent digit string a is limited to 
   * nine digits and 80 total characters, including blanks.  
   * The exponent portion (i.e. the portion up to and including x) 
   * and the period may optionally
   * be omitted.
   * Debug output starts with debug_level = 7.
   *
   * This routine has not been thoroughly tested.
   */
  int i, j, is, id; 
  double bi;
  char ai;
  int n6 = prec_words+6;
  mp_real f(0.0, 9), sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  int prec_words = mp::prec_words;
  
  if (error_no != 0) {
    if (error_no == 99)  mpabrt();
    zero(b);
    return;
  }
  
  if (debug_level >= 7) {
    int no = std::min (n, int (7.225 * debug_words) + 20);
    cerr << "MPINPC I ";
    for(i = 0; i < no; ++i) {
      cerr << a[i];
      if (i%78 == 0) cerr << endl;
    } 
  }
  
  char *ca = new char[81];
  //int nws = prec_words++;
  int i1 = 0;
  int nn = 0;
  
  //  Find the carat, period, plus or minus sign, whichever comes first.
  
  bool caretFound = false;
  
  for (i = 0; i < n; ++i) {
    ai = a[i];
    if (ai == '^') {
      caretFound = true;
      break;
    }
    if (ai == '.' || ai == '+' || ai == '-') break;
  }
  
  for (j = 0; j < 81; ++j) ca[j]='\0';
	
  if (caretFound) {
    // Make sure number preceding the caret is 10.
    int i2 = i-1;
    if (i2 > 79) {
      delete [] ca;
      mpinpcExit();
      return;
    }
		
    j = 0;
    for (i = 0; i <= i2; ++i) {
      ai = a[i];
      if (ai == ' ') continue;
      else if (!isdigit(ai)) {
        delete [] ca;
        mpinpcExit();
        return;
      }
      ca[j++] = ai;
    }

    if (ca[0]!='1' || ca[1]!='0') {
      delete [] ca;
      mpinpcExit();
      return;
    }

    // Find the x or *.
    i1 = i2 + 2; // first char after carat
    bool exit = true;
    for (i = i1; i < n; ++i) {
      ai = a[i];
      if (ai == 'x' || ai == '*') {
        exit = false;
        break;
      }
    }
    if (exit) {
      delete [] ca;
      mpinpcExit();
      return;
    }
		
    //  Convert the exponent.
    i2 = i - 1;
    int l1 = i2 - i1;
    if (l1 > 79) {
      delete [] ca;
      mpinpcExit();
      return; 
    }
    id = 0;
    is = 1;
    j = 0;
    for (i = 0; i <= l1; ++i) {
      ai = a[i+i1];
      if (ai == ' ' || ai == '+') continue;
      else if (ai == '-' && id == 0) {
        id = 1;
        is = -1;
      } else {
        if (!isdigit(ai)) {
          delete [] ca;
          mpinpcExit();
          return;
        }
        id = 1; 
        ca[j++] = ai;
      }
    } // end for i = ...

    ca[j]='\0';
    nn = atoi(ca);
    nn = is * nn;
    i1 = i2 + 2;
  } // end if (caretFound) ...

  //  Find the next nonblank character.
  bool exit = true;
  for (i = i1; i < n; ++i) {
    if (a[i] != ' ') {
      exit = false;
      break;
    }
  }
  if (exit) {
    delete [] ca;
    mpinpcExit();
    return;
  }
  
  //  Check if the nonblank character is a plus or minus sign.
  i1 = i;
  if (a[i1] == '+') {
    i1 = i1 + 1;
    is = 1;
  } else if (a[i1] == '-') {
    i1 = i1 + 1;
    is = -1;
  } else {
    is = 1;
  }

  int nb = 0;
  int ib = 0;
  id = 0;
  zero(sk2);
  f[1] = 1.;
  f[2] = 0.;
  int ip;     // position of period
  int it = 0; // iteration count
  
  int mm;
  bool cont = true;
  while (cont) {
    ip = 0;
    for (mm = 0; mm < 6; ++mm) ca[mm]='0';
		
    //  Scan for digits, looking for the period also. On the first pass we just
    //  count, so that on the second pass it will come out right.
    for (i = i1; i < n; ++i) {
      ai = a[i];
      if (ai == ' ') {
      } else if (ai == '.') {
        if (ip != 0) {
          delete [] ca;
          mpinpcExit();
          return;
        }
        ip = id; // period removed in the 1st pass, but position remembered
      } else if(ai == ',' || ai == '\t' || ai == '\r' || ai == '\n') {
	ai = ' ';
      } else if (!isdigit(ai)) {
        delete [] ca;
        mpinpcExit();
        return;
      } else {
        id++;
        ca[ib++] = ai;
      }
      if (ib == 6 || (i == (n-1) && ib != 0)) {
        if (it != 0) { // second pass
          nb++;
          ca[ib]='\0';
          
          bi = atoi(ca);
					
          mpmuld (sk2, 1e6, 0, sk0, prec_words);
          
          if (bi != 0) {
	    f[1] = 1.;
            f[3] = double(bi); // *cast*
          } else {
	    f[1] = 0.;
          }
          mpadd (sk0, f, sk2, prec_words);
          for (mm = 0; mm < 6; ++mm) ca[mm]='0';
        }
        if ( (i+1) != n ) ib = 0;
      }
    } // for i = ...
    
    if (it == 0) {
      ib = 6 - ib;
      if (ib == 6) ib = 0; // digits are multiple of 6
      it = 1;
    } else
      cont = false;
  }
  
  if (is == -1) sk2[1] = -sk2[1];
  if (ip == 0) ip = id;
  
  nn = nn + ip - id;
  f[1] = 1.;
  f[3] = 10.;
  mpnpwr(f, nn, sk0, prec_words);
  mpmul(sk2, sk0, sk1, prec_words);
  mpeq(sk1, b, prec_words);
  mproun(b);
  
  delete [] ca;
  if (debug_level >= 7) print_mpreal("MPINPC O ", b);
  return;
}

void mp_real::mpinpcExit()
{
  /**
   * A helper function for mpinpc.
   */
  if (MPKER[41] != 0) {
    cerr << "*** MPINPC: Syntax error in literal string." << endl;
    error_no = 41;
    if (MPKER[error_no] == 2)  mpabrt();
  }
}

