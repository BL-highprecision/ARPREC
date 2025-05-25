#define _CRT_SECURE_NO_DEPRECATE
#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include "small_inline.h"
#include <cstdio>

using std::cerr;
using std::endl;

/* This routine is deprecated.  
 * Use more C++ friendly mp_real::to_string or mp_real::write.  */

void mp_real::mpoutx(const mp_real& a, int la, char *b, int& n, int prec_words)
{
  /**
   * Converts the MP number A into character form in the
   * char array B.  N (an output parameter) is the length of the output.
   * In other words, B is contained in B[0]...B[N-1].  The format
   * is analogous to the Fortran exponential format (E format),
   * except that the exponent is placed first.
   * 
   * before calling MPOUTX, the arrays mpuu1 and mpuu2 must be initialized
   * by calling MPINIX.  For modest levels of precision, 
   * use MPOUTC.
   *
   * Space needed in b :   mpnbt * log_10(2) * prec_words + 30 cells
   *
   *
   * The algorithm is to multiply by a large power of ten, then
   * round down to the nearest integer, then print the integer
   * into a string (using mpoutx_help), the use the integer
   * to write the result into b.
   *
   *  or in other words 15.06 * prec_words + 30 cells.
   */
  const double al2 = 0.301029995663981195;

  if(error_no != 0) {
    if(error_no == 99) mpabrt();
    b[0] = ' ';
    b[1] = '\0';
    return;
  }
  if(debug_level >= 7) cerr << "MPOUTX I" << endl;

  int nws = prec_words;
  int ncr = 1 << (mpmcrx+1);
#if 0
  printf("mpoutx[1] prec_words %d, ncr %d\n", prec_words, ncr);
#endif
  if(prec_words < ncr || ((a[2]+1.0 >= std::abs(a[1])) && (std::abs(a[1]) < ncr))) {
    mpout(a, std::min(la, int(15.06 * std::abs(a[1]))+4), b, n, prec_words);
    return;
  }

  prec_words += 6;
  int n6 = prec_words+6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6);
  mp_real sk3(0.0, n6), sk4(0.0, n6);

  // Normalize input to an 
  // approximate integer by multiplying by a suitable power of 10.

  double t1 = a[FST_M] + mprdx * a[FST_M+1] + mprx2 * a[FST_M+2];
  double t2 = log10(t1);
  int m1;
  m1 = la - int(mpnbt*al2*a[2]) -int(t2);
  int m2 = int(((std::abs(a[1]) - a[2] - 1.0)*mpnbt));
  if(m2 > 0)
    m1 = std::min(m1, m2);
  m1+=2;

  mpdmc(10.0, 0, sk0, prec_words);
  mpnpwx(sk0, m1, sk2, prec_words);
  mpmulx(a, sk2, sk1, prec_words);
  sk1[1] = std::abs(sk1[1]);
  mpnint(sk1, sk1, prec_words);
  if((sk1[1] * al2 * mpnbt) >= la + 50) {
    cerr << "*** MPOUTX : Exponent Error : sk1 too large" << endl;
    mpabrt();
  }

  //now we have an integer, close to 10^m1 * a.
  char *integer_string = new char[la+50];
  mpoutx_help(sk1, integer_string, la+50, prec_words);
  
  int i=0, j;
  while(integer_string[i] == '0') i++;
  //i holds the index of the first nonzero digit.
  long exponent = la+50 - i - 1 - m1;
  b[0] = '1';
  b[1] = '0';
  b[2] = ' ';
  b[3] = '^';
  sprintf(b+4, "%11ld", exponent);
  b[15] = ' ';
  b[16] = 'x';
  b[17] = ' ';
  b[18] = a[1] >= 0.0 ? ' ' : '-';
  b[19] = integer_string[i++];
  b[20] = '.';
  for(j=i;j<la+50; j++) {
    b[21+j-i] = integer_string[j];
  }
  if(j-i > n_output_digits) { // got too many digits.
    while(j-i > n_output_digits) j--;
  }
  b[21+j-i] = ',';
  b[21+j-i+1] = '\0';
  n = 21+j-i+1;
  i = n-2;
  while(b[i] == '9') i--;
  if(i != n-2) {
    if(b[i] == '.') {//nines all the way to the decimal point.
      if(b[i-1] == '9') {//first digit was also nine.
	exponent++;
	sprintf(b+4, "%11ld", exponent);
	b[15] = ' ';
	b[16] = 'x';
	b[17] = ' ';
	b[18] = a[1] >= 0.0 ? ' ' : '-';
	b[i-1] = '1';
	b[i+1] = ',';
	b[i+2] = '\0';
	n = i+2;
      } else {
	b[i-1] += 1;
	b[i+1] = ',';
	b[i+2] = '\0';
	n = i+2;
      }
    } else {
      b[i] += 1;
      b[i+1] = ',';
      b[i+2] = '\0';
      n = i+2;
    }
  } else {
    i = n-2;
    while(b[i] == '0') i--;
    if(i != n-2) {
      b[i+1] = ',';
      b[i+2] = '\0';
      n = i+2;
    }
  }

  delete [] integer_string;
  prec_words = nws;
  return;  
}

void mp_real::mpoutx_help(const mp_real& a, char *b, int n, int prec_words)
{
  /**
   * converts the POSITIVE INTEGER a into base-10 and places the
   * result in b, with length n.
   *
   * leading zeros are placed in the array.
   */
  const double al2 = 0.301029995663981195;

  if(a[1] == 0.0) {
    //zero out our section.
    int i;
    for(i=0;i<n;i++) b[i] = '0';
    return;
  }

  int nws = prec_words;
  int ncr = 1 << (mpmcrx+1);

  prec_words = std::min(std::max(int(a[1])+3, int(a[2])+4), prec_words);
  //Check if actual precision level of argument is too low to 
  //justify the advanced routine.
  
  if(prec_words < ncr || ((a[2]+1.0 >= std::abs(a[1])) && (std::abs(a[1]) < ncr))) {
    char *temp_string = new char[n+200];
    int curr_output_prec;
    int i, j, n2= n;
    curr_output_prec = n_output_digits;
    n_output_digits = std::min(n_output_digits, n+10);
    //mpoutc has output of the form 10 ^ # x D.########...
    // where # represents one or more digits, and D represents
    // a single digit

    mpoutc(a, temp_string, n2, prec_words);
    if(n2 >= n+200) {
      cerr << "*** MPOUTX HELPER: non-integer or other bad input to helper" << endl;
      mpabrt();
    }
    if(n2 < n+100) {
      // add a comma.  sometimes mpout doesn't put one on.
      temp_string[n2] = ',';
      n2++;
      temp_string[n2] = '\0';
    }
    //now strip out the exponent.
    i=0;
    while(i<n+200 && temp_string[i++] != '^');
    j = i;
    while(j<n+200 && temp_string[j++] != 'x');
    if(j == n+200 || temp_string[--j] != 'x') {
      cerr << "*** MPOUTX HELPER: bad data from MPOUTC" << endl;
      mpabrt();
    }
    //now temp_string[i] is just after the ^, and
    // temp_string[j] is the x 
    temp_string[j] = '\0';
    long exponent = atol(&temp_string[i]);
    if(exponent > n) {
      cerr << "*** MPOUTX HELPER: not enough space for base-10 output in small case." << endl;
      mpabrt();
    }
    j++;
    while(j<n+200 && temp_string[j++] != '.');
    if(j == n+200 || temp_string[--j] != '.') {
      cerr << "*** MPOUTX: bad data from MPOUTC" << endl;
      mpabrt();
    }
    // now temp_string[j] is on the decimal point, just
    // after the first digit.
    // We need to decide from the exponent how many leading zeroes
    // to put down.
    int leading_zeros = n - exponent-1;
    for(i=0;i<leading_zeros;i++) b[i] = '0';
    //get first digit.
    b[i] = temp_string[j-1];
    j++; i++;
    for(;i<n;i++, j++) {    
      b[i] = temp_string[j];
      if(b[i] == ',') break;
      if(!isdigit(b[i]))
	b[i] = '0';
    }
    //sometimes mpout prints out integers incorrectly.
    if(j<n2 && 
       isdigit(temp_string[j]) && 
       temp_string[j] >= '5') { //round up.
      int i2 = i-1;
      while(i2>=0 && b[i2] == '9') i2--;
      if(i2 < 0) {
	//all nines - wierd.
	cerr << "*** MPOUTX round failed. use mpout." << endl;
	mpabrt();
      }
      b[i2]++;
      i = i2+1;
    }
    //fill in trailing zeros
    for(;i<n;i++, j++) b[i] = '0';

    //we should be done with this recurse.
    prec_words = nws;
    n_output_digits = curr_output_prec;
    delete [] temp_string;
    return;
  }

  // there are too many digits to efficiently do with
  // mpoutc.  Prepare to recurse.

  prec_words++;
  int n6 = prec_words+6;
  mp_real sk0(0.0, n6), sk1(0.0, n6), sk2(0.0, n6), sk3(0.0, n6);

  int i1, i2;
  double t1;
  // split large integer into two approximateley equal decimal sections.
  mpmdc(a, t1, i1, prec_words);
  i2 = int(log10(t1))+1;
  i2 += int(i1 * al2)+1;
  int m2 = (i2+1) / 2;
  m2++;
  if(m2 > n) {
    cerr << "*** MPOUTX HELPER: not enough space for base-10 output." << endl;
    mpabrt();
  }
  //do the split into superdigits of base 10^m2
  mpdmc(10.0, 0, sk0, prec_words);
  mpnpwx(sk0, m2, sk3, prec_words);

  //perform a mod 10^m2.  division is performed at less precision.
  int nws2 = prec_words;
  prec_words = std::min(prec_words, int(a[2]) - int(sk3[2]) + 4);
  prec_words = std::max(prec_words, 1);
  mpdivx(a, sk3, sk0, prec_words);
  prec_words = nws2;

  mpinfr(sk0, sk2, sk0, prec_words); //fractional part is garbage
  mpmulx(sk2, sk3, sk0, prec_words);
  mpsub(a, sk0, sk0, prec_words);
  //test for error in division and modulus, and fix it.
  if(sk0[1] < 0.0) {
    mpsub(sk2, mp_real(1.0), sk2, prec_words);
    mpadd(sk0, sk3, sk0, prec_words);
  } else if(sk0 >= sk3) {
    mpadd(sk2, mp_real(1.0), sk2, prec_words);
    mpsub(sk0, sk3, sk0, prec_words);
  }

  //now sk0 holds the low superdigit, and sk2 holds the high superdigit
  //(superdigits are in base 10^m2.)
  

  // Recursively convert each section.

  mpoutx_help(sk2, b, n-m2, prec_words);

  //now the other recursion
  mpoutx_help(sk0, b+(n-m2), m2, prec_words);

  prec_words = nws;
  return;
}

