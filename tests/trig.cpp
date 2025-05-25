#include <iostream>
#include <iomanip>

#include <arprec/mp_real.h>

using std::cerr;
using std::cout;
using std::endl;

bool flag_verbose = false;

bool check_result(const char *msg, const mp_real &x_comp, 
                  const mp_real &x_truth, const mp_real &eps = mp_real::_eps) {
  cout << msg << "... ";
  mp_real err = abs(x_comp - x_truth);
  bool pass = err <= eps;
  cout << (pass ? "pass" : "fail") << endl;
  if (flag_verbose || !pass) {
    cout << std::scientific << std::setprecision(32);
    cout << "  computed = " << x_comp << endl;
    cout << "  truth    = " << x_truth << endl;
    cout << "  err      = " << err << endl;
  }
  return pass;
}

void print_usage() {
  cout << "trig [-h] [-n N] [-v]" << endl;
  cout << "  Performs tests of the trigonometric functions in ARPREC library." << endl;
  cout << endl;
  cout << "  -h     Prints this usage message." << endl;
  cout << "  -n N   Use N digits of precision. [defaults to 1000]" << endl;
  cout << "  -v     Print detailed information for each test." << endl;
}

int main(int argc, char *argv[]) {
  int n_digits = 1000;
  bool pass = true;

  char *arg;
  for (int i = 1; i < argc; i++) {
    arg = argv[i];
    if (strcmp(arg, "-h") == 0) {
      print_usage();
      exit(0);
    } else if (strcmp(arg, "-n") == 0) {
      if (++i < argc)
        n_digits = atoi(argv[i]);
      else
        cerr << "A number must follow after -n." << endl;
    } else if (strcmp(arg, "-v") == 0) {
      flag_verbose = true;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
      exit(1);
    }
  }

  mp::mp_init(n_digits);

  mp_real pi = mp_real::_pi;
  mp_real sqrt2 = sqrt(mp_real(2.0));
  mp_real sqrt3 = sqrt(mp_real(3.0));
  mp_real zero = 0.0;
  mp_real one = 1.0;
  mp_real x, y;

  pass &= check_result("sin(0)", sin(zero), 0.0);
  pass &= check_result("sin(pi)", sin(pi), 0.0);
  pass &= check_result("sin(pi/2)", sin(pi / 2.0), 1.0);
  pass &= check_result("sin(pi/3)", sin(pi / 3.0), sqrt3 * 0.5);
  pass &= check_result("sin(pi/4)", sin(pi / 4.0), sqrt2 * 0.5);
  pass &= check_result("sin(pi/6)", sin(pi / 6.0), 0.5);
  pass &= check_result("sin(pi + pi/3)", sin(pi + pi / 3.0), -sqrt3 * 0.5);
  pass &= check_result("sin(pi + pi/4)", sin(pi + pi / 4.0), -sqrt2 * 0.5);
  pass &= check_result("sin(pi + pi/6)", sin(pi + pi / 6.0), -0.5);
  pass &= check_result("sin(-3*pi - pi/3)", sin(-3*pi - pi / 3.0), sqrt3 * 0.5);
  pass &= check_result("sin(-3*pi - pi/4)", sin(-3*pi - pi / 4.0), sqrt2 * 0.5);
  pass &= check_result("sin(-3*pi - pi/6)", sin(-3*pi - pi / 6.0), 0.5);

  pass &= check_result("cos(0)", cos(zero), 1.0);
  pass &= check_result("cos(pi)", cos(pi), -1.0);
  pass &= check_result("cos(pi/2)", cos(pi / 2.0), 0.0);
  pass &= check_result("cos(pi/3)", cos(pi / 3.0), 0.5);
  pass &= check_result("cos(pi/4)", cos(pi / 4.0), sqrt2 * 0.5);
  pass &= check_result("cos(pi/6)", cos(pi / 6.0), sqrt3 * 0.5);
  pass &= check_result("cos(pi - pi/3)", cos(pi - pi / 3.0), -0.5);
  pass &= check_result("cos(pi - pi/4)", cos(pi - pi / 4.0), -sqrt2 * 0.5);
  pass &= check_result("cos(pi - pi/6)", cos(pi - pi / 6.0), -sqrt3 * 0.5);
  pass &= check_result("cos(-3pi + pi/3)", cos(-3*pi + pi / 3.0), -0.5);
  pass &= check_result("cos(-3pi + pi/4)", cos(-3*pi + pi / 4.0), -sqrt2 * 0.5);
  pass &= check_result("cos(-3pi + pi/6)", cos(-3*pi + pi / 6.0), -sqrt3 * 0.5);

  pass &= check_result("tan(0)", tan(zero), 0.0);
  pass &= check_result("tan(pi)", tan(pi), 0.0);
  pass &= check_result("tan(pi/3)", tan(pi / 3.0), sqrt3);
  pass &= check_result("tan(pi/4)", tan(pi / 4.0), 1.0);
  pass &= check_result("tan(pi/6)", tan(pi / 6.0), 1.0 / sqrt3);
  pass &= check_result("tan(pi - pi/3)", tan(pi - pi / 3.0), -sqrt3);
  pass &= check_result("tan(pi - pi/4)", tan(pi - pi / 4.0), -1.0);
  pass &= check_result("tan(pi - pi/6)", tan(pi - pi / 6.0), -1.0 / sqrt3);
  pass &= check_result("tan(3pi - pi/3)", tan(3*pi - pi / 3.0), -sqrt3);
  pass &= check_result("tan(3pi - pi/4)", tan(3*pi - pi / 4.0), -1.0);
  pass &= check_result("tan(3pi - pi/6)", tan(3*pi - pi / 6.0), -1.0 / sqrt3);

  x = sin(mp_real(0.1));
  y = cos(mp_real(0.1));
  pass &= check_result("sin(0.1)^2 + cos(0.1)^2", x*x + y*y, 1.0);

  x = sin(mp_real(10.0));
  y = cos(mp_real(10.0));
  pass &= check_result("sin(100)^2 + cos(100)^2", x*x + y*y, 1.0);

  pass &= check_result("asin(0)", asin(zero), 0.0);
  pass &= check_result("asin(1)", asin(one), pi / 2.0);
  pass &= check_result("asin(sqrt(2)/2)", asin(0.5 * sqrt2), pi / 4.0);
  pass &= check_result("asin(-sqrt(3)/2)", asin(-0.5 * sqrt3), -pi / 3.0);

  pass &= check_result("acos(0)", acos(zero), pi / 2.0);
  pass &= check_result("acos(1)", acos(one), 0.0);
  pass &= check_result("acos(sqrt(2)/2)", acos(0.5 * sqrt2), pi / 4.0);
  pass &= check_result("acos(-sqrt(3)/2)", acos(-0.5 * sqrt3), pi-pi / 6.0);

  pass &= check_result("atan(0)", atan(zero), 0.0);
  pass &= check_result("atan(1)", atan(one), pi / 4.0);
  pass &= check_result("atan(sqrt(3))", atan(sqrt3), pi / 3.0);
  pass &= check_result("atan(-1/sqrt(3))", atan(-1.0 / sqrt3), -pi / 6.0);

  return pass;
}
