/*
 * tests/mp_timer.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2004
 *
 * Contains code to time basic operations.
 */

#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <arprec/mp_real.h>
#include "tictoc.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;
using std::strcmp;
using std::setprecision;
using std::setw;

// Global flags passed to the main program.
static bool flag_test_double = false;
static bool flag_test_mp = false;
bool flag_verbose = false;
int nr_digits = 1000;


template <class T>
class TestSuite {
public:
  void test1();
  void test2();
  void test3();
  void test4();
  void test5();
  void test6();
  void testall();
  T pi();
};

template <>
mp_real TestSuite<mp_real>::pi() { return mp_real::_pi; }

template <>
double TestSuite<double>::pi() { return 3.141592653589793116; }

void print_timing(double nops, double t) {
  double mops = 1.0e-3 * nops / t;
  if (flag_verbose) {
    cout << std::scientific << setprecision(3);
    cout << "  " << nops << " operations in " << t << " seconds." << endl;
    cout << "  " << 1.0/mops << " ms/op  " << mops << " kop/s." << endl;
  } else {
    cout << std::fixed;
    cout << setprecision(6) << setw(10) << 1.0 / mops << " ms";
    cout << setprecision(4) << setw(10) << mops << " kop/s" << endl;
  }
}

template <class T>
void TestSuite<T>::test1() {
  if (flag_verbose)
    cout << endl << "Timing addition..." << endl;
  else
    cout << "  add: ";

  int n = 20000;
  tictoc_t tv;
  double t;

  T a1 = pi() * 1.0e-100;
  T a2 = pi();
  T b1 = 1.0, b2 = 0.0;

  tic(&tv);
  for (int i = 0; i < n; i++) {
    b1 += a1;
    b2 += a2;
  }
  t = toc(&tv);

  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "(b1-1)/pi-1 = " << (b1-1.0)/(pi()*1e-100) - 1.0 << endl;
    cout << "b2/pi-1     = " << b2/pi() - 1.0 << endl;
  }
  print_timing(2*n, t);
}

template <class T>
void TestSuite<T>::test2() {
  if (flag_verbose)
    cout << endl << "Timing multiplication ..." << endl;
  else
    cout << "  mul: ";

  int n = 2000;
  tictoc_t tv;
  double t;

  T a1 = 1.0 + 0.125 * pi() / n;
  T a2 = 1.0 - pi() / n;
  T b1 = 1.0, b2 = 1.0;

  tic(&tv);
  for (int i = 0; i < n; i++) {
    b1 *= a1;
    b2 *= a2;
  }
  t = toc(&tv);

  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "b1-exp(pi/8) = " << b1-exp(0.125*pi()) << endl;
    cout << "b2-exp(-pi)  = " << b2-exp(-pi()) << endl;
  }
  print_timing(2*n, t);
}

template <class T>
void TestSuite<T>::test3() {
  if (flag_verbose)
    cout << endl << "Timing division ..." << endl;
  else
    cout << "  div: ";

  int n = 2000;
  tictoc_t tv;
  double t;

  T a1 = 1.0 + pi() / n;
  T a2 = 1.0 - 0.125 * pi() / n;
  T b1 = 1.0, b2 = 1.0;

  tic(&tv);
  for (int i = 0; i < n; i++) {
    b1 /= a1;
    b2 /= a2;
  }
  t = toc(&tv);

  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "b1-exp(-pi)  = " << b1-exp(-pi()) << endl;
    cout << "b2-exp(pi/8) = " << b2-exp(0.125*pi()) << endl;
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test4() {
  if (flag_verbose)
    cout << endl << "Timing square root ..." << endl;
  else
    cout << " sqrt: ";

  int n = 500;
  tictoc_t tv;
  double t;

  T a = 0.0;

  tic(&tv);
  for (int i = 0; i < n; i++) {
    a = sqrt(2.0 + a);
  }
  a = sqrt(2.0 - a);
  t = toc(&tv);
  for (int i = 0; i <= n; i++) {
    a *= 2.0;
  }
  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "a-pi = " << a-pi() << endl;
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test5() {
  if (flag_verbose)
    cout << endl << "Timing sin ..." << endl;
  else 
    cout << "  sin: ";

  int n = 100;
  tictoc_t tv;
  double t;

  T a = 0.0;
  T c = 1.7 * T(1.0) / double(n);
  T d = 2.45 * pi() / double(n + 3);

  tic(&tv);
  for (int i = 0; i < n; i++) {
    a = a + sin(c);
    c += d;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "a = " << a << endl;
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::test6() {
  if (flag_verbose)
    cout << endl << "Timing log ..." << endl;
  else
    cout << "  log: ";

  int n = 100;
  tictoc_t tv;
  double t;

  T a = 0.0;
  T c = exp(T(-50.1));
  T d = exp(T(100.2) / double(n));

  tic(&tv);
  for (int i = 0; i < n; i++) {
    a = a + log(c);
    c *= d;
  }
  t = toc(&tv);
  if (flag_verbose) {
    cout << std::scientific << setprecision(10);
    cout << "a = " << a << endl;
  }

  print_timing(n, t);
}

template <class T>
void TestSuite<T>::testall() {
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
}

void print_usage() {
  cout << "mp_test [-h] [-double] [-mp] [-all]" << endl;
  cout << "  Performs timing tests of the arprec library." << endl;
  cout << endl;
  cout << "  -h -help  Prints this usage message." << endl;
  cout << "  -double   Time arithmetic of double." << endl;
  cout << "  -mp       Time arithmetic of mp_real arithmetic." << endl;
  cout << "  -all      Time both double arithmetic and mp_real arithmetic." << endl;
  cout << "  -n N      Use N digits of precision.  (default is 1000)." << endl;
  cout << "  -v        Verbose output." << endl;
}

int main(int argc, char *argv[]) {

  /* Parse the arguments. */
  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0 || strcmp(arg, "-help") == 0) {
      print_usage();
      std::exit(0);
    } else if (strcmp(arg, "-double") == 0) {
      flag_test_double = true;
    } else if (strcmp(arg, "-mp") == 0) {
      flag_test_mp = true;
    } else if (strcmp(arg, "-n") == 0) {
      if (++i < argc)
        nr_digits = atoi(argv[i]);
      else {
        cerr << "A number must follow after -n." << endl;
        exit(-1);
      }
    } else if (strcmp(arg, "-all") == 0) {
      flag_test_double = flag_test_mp = true;
    } else if (strcmp(arg, "-v") == 0) {
      flag_verbose = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  mp::mp_init(nr_digits);

  /* If no flag, test both double-double and quad-double. */
  if (!flag_test_double && !flag_test_mp) {
    flag_test_mp = true;
  }

  if (flag_test_double) {
    TestSuite<double> test;

    cout << endl;
    cout << "Timing double" << endl;
    test.testall();
  }

  if (flag_test_mp) {
    TestSuite<mp_real> test;

    cout << endl;
    cout << "Timing mp_real (" << nr_digits << " digits)" << endl;
    test.testall();
  }

  mp::mp_finalize();
  return 0;
}

