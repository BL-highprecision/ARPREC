#include <iostream>
#include <iomanip>
#include <arprec/mp_real.h>

using std::cout;
using std::cerr;
using std::endl;

// Global flags passed to the main program.
bool flag_verbose = false;
int nr_digits = 1000;

mp_real polyeval(const mp_real *c, int n, const mp_real &x) {
  /* Just use Horner's method of polynomial evaluation. */
  mp_real r = c[n];

  for (int i = n-1; i >= 0; i--) {
    r *= x;
    r += c[i];
  }

  return r;
}

mp_real polyroot(const mp_real *c, int n, 
                 const mp_real &x0, const mp_real &eps) {
  mp_real x = x0;
  mp_real f;
  mp_real *d = new mp_real[n];
  bool conv = false;
  int i;

  /* Compute the coefficients of the derivatives. */
  for (i = 0; i < n; i++) {
    d[i] = c[i+1] * (i+1);
  }

  /* Newton iteration. */
  for (i = 0; i < 100; i++) {
    f = polyeval(c, n, x);

    if (abs(f) < eps) {
      conv = true;
      break;
    }
    x -= (f / polyeval(d, n-1, x));
  }
  delete [] d;

  if (!conv) {
    cerr << "ERROR (dd_real::polyroot): Failed to converge." << endl;
    exit(-1);
    return 0.0;
  }

  return x;

}

bool print_result(bool result) {
  if (result)
    cout << "Test passed." << endl;
  else
    cout << "Test FAILED." << endl;
  return result;
}

template <class T>
class TestSuite {
  mp_real _one;
  mp_real _eps;
  mp_real _pi;
  mp_real _e;
  mp_real _log2;
public:
  TestSuite(int nr_digits) {
    _one = 1.0;
    _eps = pow(mp_real(10.0), -nr_digits+10);
    mp_real::mppi(_pi);
    _e = exp(_one);
    _log2 = log(mp_real(2.0));
  }
  bool test1();
  bool test2();
  bool test3();
  bool test4();
  bool test5();
  bool test6();
  bool testall();
};


/* Test 1.   Polynomial Evaluation / Polynomial Solving */
template <class T>
bool TestSuite<T>::test1() {
  cout << endl;
  cout << "Test 1.  (Polynomial)." << endl;

  static const int n = 8;
  T *c = new T[n];
  T x, y;

  for (int i = 0; i < n; i++)
    c[i] = i+1;

  x = polyroot(c, n-1, 0.0, _eps);
  y = polyeval(c, n-1, x);

  if (flag_verbose) {
    cout << "Root Found:  x  = " << x << endl;
    cout << "           p(x) = " << y << endl;
  }

  delete [] c;
  return (y < 4.0*_eps);
}

/* Test 2.  Machin's Formula for Pi. */
template <class T>
bool TestSuite<T>::test2() {

  cout << endl;
  cout << "Test 2.  (Machin's Formula for Pi)." << endl;
  
  /* Use the Machin's arctangent formula:

       pi / 4  =  4 arctan(1/5) - arctan(1/239)

     The arctangent is computed based on the Taylor series expansion

       arctan(x) = x - x^3 / 3 + x^5 / 5 - x^7 / 7 + ...
  */

  T s1, s2, t, r;
  int k;
  int sign;
  mp_real d;
  mp_real err;

  /* Compute arctan(1/5) */
  d = 1.0;
  t = T(1.0) / 5.0;
  r = sqr(t);
  s1 = 0.0;
  k = 0;

  sign = 1;
  while (t > _eps) {
    k++;
    if (sign < 0)
      s1 -= (t / d);
    else
      s1 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    cout << k << " Iterations" << endl;

  /* Compute arctan(1/239) */
  d = 1.0;
  t = T(1.0) / 239.0;
  r = sqr(t);
  s2 = 0.0;
  k = 0;

  sign = 1;
  while (t > _eps) {
    k++;
    if (sign < 0)
      s2 -= (t / d);
    else
      s2 += (t / d);

    d += 2.0;
    t *= r;
    sign = -sign;
  }

  if (flag_verbose)
    cout << k << " Iterations" << endl;

  T p = 4.0 * s1 - s2;

  p *= 4.0;
  err = abs(p - _pi);

  if (flag_verbose) {
    cout << "   pi = " << p << endl;
    cout << "  _pi = " << _pi << endl;
    cout << "error = " << err << endl;
  }

  return (err < 4.0*_eps);
}

/* Test 3.  Salamin-Brent Quadratic Formula for Pi. */
template <class T>
bool TestSuite<T>::test3() {
  cout << endl;
  cout << "Test 3.  (Salamin-Brent Quadratic Formula for Pi)." << endl;

  T a, b, s, p;
  T a_new, b_new, s_new;
  mp_real m;
  mp_real err;

  a = 1.0;
  b = sqrt(T(0.5));
  s = 0.5;
  m = 1.0;

  p = 2.0 * sqr(a) / s;
  if (flag_verbose)
    cout << "Iteration 0: " << p << endl;
  for (int i = 1; i <= 20; i++) {
    m *= 2.0;
    a_new = 0.5 * (a + b);
    b_new = sqrt(a * b);
    s_new = s - m * (sqr(a_new) - sqr(b_new));
    a = a_new;
    b = b_new;
    s = s_new;
    p = 2.0 * sqr(a) / s;
    if (flag_verbose)
      cout << "Iteration " << i << ": " << p << endl;
  }

  err = abs(p - _pi);

  if (flag_verbose) {
    cout << "        _pi: " << _pi << endl;
    cout << "      error: " << err << " = " << err / _eps << " eps" << endl;
  }

  // for some reason, this test gives relatively large error compared
  // to other tests.  May need to be looked at more closely.
  return (err < 1024.0*_eps);
}

/* Test 4.  Borwein Quartic Formula for Pi. */
template <class T>
bool TestSuite<T>::test4() {
  cout << endl;
  cout << "Test 4.  (Borwein Quartic Formula for Pi)." << endl;

  T a, y, p, r;
  mp_real m;
  mp_real err;

  a = 6.0 - 4.0 * sqrt(T(2.0));
  y = sqrt(T(2.0)) - 1.0;
  m = 2.0;

  p = 1.0 / a;
  if (flag_verbose)
    cout << "Iteration 0: " << p << endl;

  for (int i = 1; i <= 9; i++) {
    m *= 4.0;
    mp_real::mpnrt(1.0 - sqr(sqr(y)), 4, r, mp::prec_words);
    y = (1.0 - r) / (1.0 + r);
    a = a * sqr(sqr(1.0 + y)) - m * y * (1.0 + y + sqr(y));
    
    p = 1.0 / a;
    if (flag_verbose)
      cout << "Iteration " << i << ": " << p << endl;
  }

  err = abs(p - _pi);
  if (flag_verbose) {
    cout << "        _pi: " << _pi << endl;
    cout << "      error: " << err << " = " << err / _eps << " eps" << endl;
  }  

  return (err < 64.0*_eps);
}

/* Test 5.  Taylor Series Formula for E. */
template <class T>
bool TestSuite<T>::test5() {

  cout << endl;
  cout << "Test 5.  (Taylor Series Formula for E)." << endl;

  /* Use Taylor series

       e = 1 + 1 + 1/2! + 1/3! + 1/4! + ...

     To compute e.
  */

  T s = 2.0, t = 1.0;
  mp_real n = 1.0;
  mp_real delta;
  int i = 0;

  while (t > _eps) {
    i++;
    n += 1.0;
    t /= n;
    s += t;
  }

  delta = abs(s - _e);

  if (flag_verbose) {
    cout << "    e = " << s << endl;
    cout << "   _e = " << _e << endl;
    cout << "error = " << delta << endl;
    cout << i << " iterations." << endl;
  }

  return (delta < 4.0*_eps);
}

/* Test 6.  Taylor Series Formula for log 2.*/
template <class T>
bool TestSuite<T>::test6() {
  cout << endl;
  cout << "Test 6.  (Taylor Series Formula for Log 2)." << endl;

  /* Use the Taylor series

      -log(1-x) = x + x^2/2 + x^3/3 + x^4/4 + ...

     with x = 1/2 to get  log(1/2) = -log 2.
  */

  T s = 0.5;
  T t = 0.5;
  mp_real delta;
  double n = 1.0;
  double i = 0;

  while (abs(t) > _eps) {
    i++;
    n += 1.0;
    t *= 0.5;
    s += (t/n);
  }

  delta = abs(s - _log2);

  if (flag_verbose) {
    cout << " log2 = " << s << endl;
    cout << "_log2 = " << _log2 << endl;
    cout << "error = " << delta << endl;
    cout << i << " iterations." << endl;
  }

  return (delta < 4.0*_eps);
}

template <class T>
bool TestSuite<T>::testall() {
  bool pass = true;
  pass &= print_result(test1());
  pass &= print_result(test2());
  pass &= print_result(test3());
  pass &= print_result(test4());
  pass &= print_result(test5());
  pass &= print_result(test6());
  return pass;
}

void print_usage() {
  cout << "mp_test [-h] [-v]" << endl;
  cout << "  Performs miscellaneous tests of the ARPREC library," << endl;
  cout << "  such as polynomial root finding, computation of pi, etc." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -v" << endl;
  cout << "  -verbose  Print detailed information for each test." << endl;
  cout << "  -n N  Use N digits of precision. (default is 1000). " << endl;
}

int main(int argc, char *argv[]) {
  
  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      exit(0);
    } else if (strcmp(arg, "-n") == 0) {
      if (++i < argc)
        nr_digits = atoi(argv[i]);
      else {
        cerr << "A number must follow after -n." << endl;
        nr_digits = 100;
      }
    } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "-verbose") == 0) {
      flag_verbose = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }


  mp::mp_init(nr_digits+5);

  TestSuite<mp_real> mp_test(nr_digits);
  cout << std::setprecision(nr_digits) << endl;
  cout << "Testing mp_real ..." << endl;
  cout << "Using " << nr_digits << " digits precision." << endl;
  if (flag_verbose)
    cout << "sizeof(mp_real) = " << sizeof(mp_real) << endl;

  bool pass = mp_test.testall();

  mp::mp_finalize();
  return (pass ? 0 : 1);
}
