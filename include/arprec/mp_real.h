/*
 * include/arprec/mp_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef ARPREC_MPREAL_H
#define ARPREC_MPREAL_H

#include <iostream>
#include <string>

#include <arprec/arprec_config.h>

#ifdef ARPREC_DEBUG
#ifdef ARPREC_VALGRIND
#include <valgrind/memcheck.h>
#else
#define VALGRIND_MAKE_MEM_UNDEFINED(x, y)
#endif
#endif

#include <arprec/mp.h>
#include <arprec/mp_real_temp.h>

#if (ARPREC_QD)
#include <qd/dd_real.h>
#include <qd/qd_real.h>
#endif

/*
 *  This class represents MP real numbers.
 */
struct ARPREC_API mp_real: public mp {
  double *mpr;
  bool    alloc;

  static mp_real _pi;
  static mp_real _log2;
  static mp_real _log10;
  static mp_real _eps;

  //high precision helper functions.
  static void mplconv(int iq, int n, int nsq, double *a, 
		      double *b, double *c);
  static void mpfftrc(int is, int m, int n, int nsq, double *x, double *y);
  static void mpfft1(int is, int m, int m1, int m2, 
		     double *x, double *y);
  static void mpfft1_radix_two(int is, int m, int n1, int n2, 
		     double *x, double *y);
  static void mpfft3(int is, int l, int ns, int m, int n, 
			    double *x, double *y);
  static void mpfft3_radix2(int is, int l, int ns, int m, int n, 
			    double *x, double *y);
  static void mpfft2(int is, int ns, int m, int n, double *x, double *y,
		     double* &out);
  static void mpfftcr(int is, int m, int n, int nsq, 
		      double *x, double *y);
  static void mpoutx_help(const mp_real& a, char *b, int n, int prec_words);

  bool construct(const std::string &expn_str, 
                 const std::string &digit_str1, const std::string &digit_str2);

  //basic procedures
  static void zero(mp_real &a) { a.mpr[1] = a.mpr[2] = 0.; }
  static void mpadd(const mp_real &a, const mp_real &b, mp_real& c, int prec_words);
  static void mpeq(const mp_real& a, mp_real& b, int prec_words);
  static void mpdiv(const mp_real& a, const mp_real& b, mp_real& c, int prec_words);
  static void mpdivd(const mp_real& a, double b, int n, mp_real& c, int prec_words);
  static void mpmul(const mp_real& a, const mp_real&b, mp_real& c, int prec_words);
  static void mpmuld(const mp_real& a, double b, int n, mp_real& c, int prec_words);
  static void mpmulacc(const mp_real &a, const mp_real &b, mp_real &c);
  static void mpnorm(double d[], mp_real &a, int prec_words);
  static void mpnpwr(const mp_real& a, int n, mp_real& b, int prec_words);
  static void mpdexc(const char a[], int l, mp_real& b);
  static void mpinpc(const char a[], int n, mp_real& b);
  static void mpinpcExit();
  static void mpout(const mp_real& a, int la, char* cs, int& l, int prec_words);
  static void mpoutc(const mp_real& a, char* b, int& n, int prec_words);
  static void mproun(mp_real &a);
  static void mpsub(const mp_real &a, const mp_real &b, mp_real& c, int prec_words);
  static int mpcpr(const mp_real& a, const mp_real& b);
  static void mpexp(const mp_real& a, const mp_real& al2, mp_real& b, int prec_words);
  static void mplog(const mp_real& a, const mp_real& al2, mp_real& b,
		    int prec_words, int nit = 3);
  static void mpnint(const mp_real& a, mp_real &b, int prec_words);
  static void mpcssh(const mp_real& a, const mp_real& al2, mp_real& x,
		     mp_real& y);
  static void mpcssn(const mp_real& a, const mp_real& pi, mp_real &x,
		     mp_real &y, int prec_words);
  static void mppi(mp_real& pi);
  static void mpang(const mp_real& x, const mp_real& y, 
		      const mp_real& pi, mp_real& a);
  static void mpsqrt(const mp_real& a, mp_real& b, int prec_words);

  // high precision routines
  static void mppix(mp_real& pi);
  static void mpmulx(const mp_real& a, const mp_real& b, mp_real& c, int prec_words);
  static void mpdivx(const mp_real& a, const mp_real& b, mp_real& c, int prec_words);
  static void mpsqx(const mp_real& a, mp_real& b, int prec_words);
  static void mpsqrtx(const mp_real& a, mp_real& b, int prec_words);
  static void mpexpx(const mp_real& a, const mp_real& pi, 
		     const mp_real& al2, mp_real& b);
  static void mplogx(const mp_real& a, const mp_real& pi, const mp_real& al2, 
		     mp_real& b, int prec_words);
  static void mpagmx(mp_real &a, mp_real& b, int prec_words);
  static void mpnpwx(const mp_real& a, int n, mp_real& b, int prec_words);
  static void mpcshx(const mp_real& a, const mp_real& pi, 
		     const mp_real& al2, mp_real& x, mp_real& y);
  static void mpinix(int n); 
  static void mpoutx(const mp_real& a, int la, char *b, int& n, int prec_words);

  //basic procedures that users should have access to. 
  static void mprand(mp_real& a);
  static void mpinfr(const mp_real& a, mp_real& b, mp_real& c, 
		     int prec_words, int want_fraction=1);
  static void mpdmc(double a, int n, mp_real& b, int prec_words);
  static void mpmdc(const mp_real& a, double& d, int& n, int prec_words);
  static void mpsort(int n, mp_real *a, int *ip);
  static void mpangx(const mp_real& x, const mp_real& y, 
		     const mp_real& pi, mp_real& a);
  static void mpcssx(const mp_real& a, const mp_real& pi,
		     mp_real& x, mp_real& y);
  static void mpnrt(const mp_real& a, int n, mp_real& b, int prec_words);
  static void mpnrtx(const mp_real& a, int n, mp_real& b);
  static void mpdotd(int n, int isa, const mp_real a[],
                     int isb, const double db[], mp_real& c);

  int to_digits(char *digits, int &expn, int n) const;
  bool write(std::ostream &s, int precision = mp::n_output_digits, 
      int width = 0, int expn_width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      int n_cols = 0, bool showpos = false, bool uppercase = false, 
      char fill = ' ') const;
  bool read(std::istream &s);
  bool read(const std::string &s);
  std::string to_string(int precision = mp::n_output_digits, 
      int width = 0, int expn_width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;

  // low-level binary input/output.
  void write_binary(std::ostream &out);
  void read_binary(std::istream &in);

  // Subscript operator. For non-const objects.
  double& operator[](int index) { return mpr[index]; }

  // Subscript operator. For const objects.
  double operator[](int index) const { return mpr[index]; }

  // Assignment operator.
  mp_real& operator=(const int&);
  mp_real& operator=(const double&);
  mp_real& operator=(const mp_real&);
  mp_real& operator=(const mp_real_temp&);
  mp_real& operator=(const char*);
  mp_real& operator+=(const mp_real&);
  mp_real& operator-=(const mp_real&); 
  mp_real& operator*=(const mp_real&);
  mp_real& operator*=(double);
  mp_real& operator/=(const mp_real&);
  mp_real& operator/=(double);


//  operator double() const;
#if 0
  mp_real& operator=(const mp_complex&);
  mp_real& operator=(const mp_int&);
#endif  

  // destructor
  ~mp_real() { if(mpr && alloc) delete [] mpr; }

  // memory allocation routines
  void allocate(unsigned int sz);
  void reallocate(unsigned int sz);

  // constructors
  mp_real();
  mp_real(double a, int sz = n_words);

  mp_real(const mp_real& x);
  mp_real(const mp_real_temp& x);
  mp_real(const char *s, unsigned int sz = n_words);
  mp_real(const std::string &s, unsigned int sz = n_words);

#ifdef ARPREC_QD
  mp_real(const dd_real &a, unsigned int sz = n_words);
  mp_real(const qd_real &a, unsigned int sz = n_words);
#endif

  explicit mp_real(const double *a);
  explicit mp_real(double *a);

  /* Allocates an array of n mp_reals.  The raw data array will be
     in a contiguous block, making it easier for block transfer operations.
     The array created must be freed using free_array(mp_real *).
   */
  static mp_real *alloc_array(int n);
  static void free_array(mp_real *m);

  mp_real_temp toTempAndDestroy() { 
    mp_real_temp ret(this->mpr); mpr = 0; return ret;
  }

  // Helper function for debugging.
  static void print_mpreal(const char* str, const mp_real& mpr, 
      std::ostream &os = std::cerr);
  void print_short(const char *name) const;

};

template <class Ch, class Tr>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& s, const mp_real& a) {
  bool showpos = (s.flags() & std::ios_base::showpos) != 0;
  bool uppercase =  (s.flags() & std::ios_base::uppercase) != 0;
  std::ios_base::fmtflags fmt = s.flags();

  a.write(s, s.precision(), s.width(), 0, fmt, 0, showpos, uppercase, s.fill());
  return s;
}

template <class Ch, class Tr>
std::basic_istream<Ch, Tr>& operator>>(std::basic_istream<Ch, Tr>& s, mp_real& a) {
  a.read(s);
  return s;
}

std::string to_string(mp_real &x, int precision = mp::n_output_digits, 
    int width = 0, int expn_width = 0, 
    std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
    bool showpos = false, bool uppercase = false, char fill = ' ');

ARPREC_API double dble(const mp_real& a);

#if (ARPREC_QD)
ARPREC_API dd_real to_dd_real(const mp_real &a);
ARPREC_API qd_real to_qd_real(const mp_real &a);
#endif

ARPREC_API mp_real_temp operator-(const mp_real &a, const mp_real &b);
ARPREC_API mp_real_temp operator+(const mp_real &a, const mp_real &b);
ARPREC_API mp_real_temp operator*(const mp_real &a, const mp_real &b);
ARPREC_API mp_real_temp operator*(const mp_real &a, const double b);
ARPREC_API mp_real_temp operator*(const double b, const mp_real& a);
ARPREC_API mp_real_temp operator/(const mp_real &a, const mp_real &b);
ARPREC_API mp_real_temp operator/(const mp_real &a, const double b);
ARPREC_API mp_real_temp operator/(const double b, const mp_real &a);

ARPREC_API mp_real_temp operator-(const mp_real &a);

ARPREC_API bool operator>(const mp_real& a, const mp_real& b);
ARPREC_API bool operator>=(const mp_real& a, const mp_real& b);
ARPREC_API bool operator<(const mp_real& a, const mp_real& b);
ARPREC_API bool operator<=(const mp_real& a, const mp_real& b);
ARPREC_API bool operator==(const mp_real& a, const mp_real& b);
ARPREC_API bool operator!=(const mp_real& a, const mp_real& b);

ARPREC_API bool operator>(double a, const mp_real& b);
ARPREC_API bool operator>=(double a, const mp_real& b);
ARPREC_API bool operator<(double a, const mp_real& b);
ARPREC_API bool operator<=(double a, const mp_real& b);
ARPREC_API bool operator==(double a, const mp_real& b);
ARPREC_API bool operator!=(double a, const mp_real& b);

ARPREC_API bool operator>(const mp_real& a, double b);
ARPREC_API bool operator>=(const mp_real& a, double b);
ARPREC_API bool operator<(const mp_real& a, double b);
ARPREC_API bool operator<=(const mp_real& a, double b);
ARPREC_API bool operator==(const mp_real& a, double b);
ARPREC_API bool operator!=(const mp_real& a, double b);

ARPREC_API mp_real_temp pow(const mp_real& a, int n);
ARPREC_API mp_real_temp pow(const mp_real& a, const mp_real& b);
ARPREC_API mp_real_temp pow(const mp_real& a, double b);
ARPREC_API mp_real_temp abs(const mp_real& a);
ARPREC_API mp_real_temp sqrt(const mp_real& a);
ARPREC_API mp_real_temp aint(const mp_real& a);
ARPREC_API mp_real_temp anint(const mp_real& a);
ARPREC_API mp_real_temp exp(const mp_real& a);
ARPREC_API mp_real_temp log(const mp_real& a);
ARPREC_API mp_real_temp log10(const mp_real& a);
ARPREC_API mp_real_temp sin(const mp_real& a);
ARPREC_API mp_real_temp cos(const mp_real& a);
ARPREC_API mp_real_temp tan(const mp_real& a);
ARPREC_API mp_real_temp sinh(const mp_real& a);
ARPREC_API mp_real_temp cosh(const mp_real& a);
ARPREC_API mp_real_temp tanh(const mp_real& a);
ARPREC_API mp_real_temp atan(const mp_real& a);
ARPREC_API mp_real_temp atan2(const mp_real& y, const mp_real& x);
ARPREC_API mp_real_temp asin(const mp_real& a);
ARPREC_API mp_real_temp acos(const mp_real& a);
ARPREC_API mp_real_temp sqr(const mp_real& a);/*square of a*/
ARPREC_API mp_real_temp fmod(const mp_real& a, const mp_real& b);
ARPREC_API mp_real_temp mp_rand();

ARPREC_API void mpcsshf(const mp_real& a, mp_real& cosh, mp_real& sinh);
ARPREC_API void mpcssnf(const mp_real& a, mp_real& cosine, mp_real& sine);

ARPREC_API mp_real erfc(const mp_real&  arg);
ARPREC_API mp_real_temp erf(const mp_real& arg);
ARPREC_API mp_real_temp bessel(const mp_real& t);
ARPREC_API mp_real_temp besselexp(const mp_real& t);
ARPREC_API mp_real_temp gamma(const mp_real& t);

double mp_two_prod_positive(double a, double b, double &err);
double mp_two_prod(double a, double b, double &err);

template <class T>
inline T sign(T a, T b) { return (b>=0 ? std::abs(a) : -std::abs(a)); }

#if (ARPREC_INLINE)
#include <arprec/mp_inline.h>
#else
double quick_two_sum(double a, double b, double &err);
double quick_two_diff(double a, double b, double &err);
double two_sum(double a, double b, double &err);
double two_diff(double a, double b, double &err);
void split(double a, double &hi, double &lo);
double two_sqr(double a, double &err);
double nint(double d);
double SLOPPY_ANINT_POSITIVE(double a); 
double FLOOR_POSITIVE(double a);
double CEIL_POSITIVE(double a);
double AINT(double a);
double POSITIVE_AINT(double a);
double ANINT(double a); 
double SLOPPY_ANINT(double a);
double aint(double a);
double anint(double a);
void sincosh(double t, double &sinh_t, double &cosh_t);
double sqr(double t);
void dd_add_d(double a[], double b, double c[]);
void dd_add_dd(double a[], double b[], double c[]);
double mp_two_prod(double a, double b, double &err);
double mp_two_prod_positive(double a, double b, double &err);

mp_real_temp operator+(const mp_real &a, const mp_real &b);
mp_real_temp operator-(const mp_real &a, const mp_real &b);
mp_real_temp operator*(const mp_real &a, const mp_real &b);
mp_real_temp operator*(const mp_real &a, const double b);
mp_real_temp operator*(const double b, const mp_real &a);
mp_real_temp operator/(const mp_real &a, const mp_real &b);
mp_real_temp operator/(const mp_real &a, const double b);
mp_real_temp operator/(const double b, const mp_real &a);
bool operator>(const mp_real& a, const mp_real& b);
bool operator>=(const mp_real& a, const mp_real& b);
bool operator<(const mp_real& a, const mp_real& b);
bool operator<=(const mp_real& a, const mp_real& b);
bool operator==(const mp_real& a, const mp_real& b);
bool operator!=(const mp_real& a, const mp_real& b);
mp_real_temp abs(const mp_real& a);
mp_real_temp acos(const mp_real& a);
mp_real_temp aint(const mp_real& a);
mp_real_temp anint(const mp_real& a);
mp_real_temp asin(const mp_real& a);
mp_real_temp atan(const mp_real& a);
mp_real_temp atan2(const mp_real& y, const mp_real& x);
mp_real_temp cos(const mp_real& a);
mp_real_temp cosh(const mp_real& a);
mp_real_temp exp(const mp_real& a);
mp_real_temp fmod(const mp_real& a, const mp_real &b);
mp_real_temp log(const mp_real& a);
mp_real_temp log10(const mp_real& a);
mp_real_temp mp_rand();
void mpcsshf(const mp_real& a, mp_real& b, mp_real& c);
void mpcssnf(const mp_real& a, mp_real& b, mp_real& c);
mp_real_temp pow(const mp_real& a, int n);
mp_real_temp pow(const mp_real& a, const mp_real& b);
mp_real_temp pow(const mp_real& a, double b);
mp_real_temp sin(const mp_real& a);
mp_real_temp sinh(const mp_real& a);
mp_real_temp sqrt(const mp_real& a);
mp_real_temp sqr(const mp_real& a);
mp_real_temp tan(const mp_real& a);
mp_real_temp tanh(const mp_real& a);
mp_real_temp operator+(const mp_real &a, int b);
mp_real_temp operator+(int b, const mp_real &a);
mp_real operator+=(mp_real &a, int b);
mp_real_temp operator-(int b, const mp_real &a);
mp_real_temp operator-(const mp_real &a, int b);
mp_real operator-=(mp_real &a, int b);
mp_real_temp operator*(const mp_real &a, int b);
mp_real operator*=(mp_real &a, int b);
mp_real_temp operator*(int b, const mp_real &a);
mp_real_temp operator/(const mp_real &a, int b);
mp_real operator/=(mp_real &a, int b);
mp_real_temp operator/(int b, const mp_real &a);
bool operator>(int b, const mp_real &a);
bool operator>(const mp_real &a, int b);
bool operator<(int b, const mp_real &a);
bool operator<(const mp_real &a, int b);
bool operator<=(int b, const mp_real &a);
bool operator<=(const mp_real &a, int b);
bool operator>=(int b, const mp_real &a);
bool operator>=(const mp_real &a, int b);
bool operator==(int b, const mp_real &a);
bool operator==(const mp_real &a, int b);
bool operator!=(int b, const mp_real &a);
bool operator!=(const mp_real &a, int b);
mp_real_temp operator+(const mp_real &a, double b);
mp_real_temp operator+(double b, const mp_real &a);
mp_real operator+=(mp_real &a, double b);
mp_real_temp operator-(double b, const mp_real &a);
mp_real_temp operator-(const mp_real &a, double b);
mp_real operator-=(mp_real &a, double b);
bool operator>(double b, const mp_real &a);
bool operator>(const mp_real &a, double b);
bool operator<(double b, const mp_real &a);
bool operator<(const mp_real &a, double b);
bool operator<=(double b, const mp_real &a);
bool operator<=(const mp_real &a, double b);
bool operator>=(double b, const mp_real &a);
bool operator>=(const mp_real &a, double b);
bool operator==(double b, const mp_real &a);
bool operator==(const mp_real &a, double b);
bool operator!=(double b, const mp_real &a);
bool operator!=(const mp_real &a, double b);
#endif

#endif
/* ARPREC_MPREAL_H */
