#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <arprec/mp_real.h>

using std::cout;
using std::endl;
using std::cerr;
using std::string;

static const char tmp_fname[] = ",,test";
static bool all_pass = true;
static int verbose = 0;

bool test_binary_io() {
  cout << "Testing binary I/O..." << endl;

  mp_real a, b, d;
  std::ofstream out(tmp_fname, std::ios_base::out | std::ios_base::binary);
  a = mp_real(1.0) / mp_real(101.0);
  a.write_binary(out);
  out.close();
  if (verbose > 0) cout << "  a written: " << a << endl;

  std::ifstream in(tmp_fname, std::ios_base::in | std::ios_base::binary);
  b.read_binary(in);
  in.close();

  if (verbose > 0) {
    cout << "  b    read: " << b << endl;
    if (a == b)
      cout << "  a equals b." << endl;
    else
      cout << "  a does not equal b." << endl;
    cout << "  a - b = " << (a-b) << endl;
  }

  std::remove(tmp_fname);
  return (a == b);
}

void process_to_digits(mp_real x, int n, 
    const char *true_digits, int true_expn, int true_len, bool &pass) {
  char buf[256];
  int expn, len;
  bool p;

  len = x.to_digits(buf, expn, n);
  p = (std::strncmp(buf, true_digits, 255) == 0) && (expn == true_expn) && (len == true_len);
  if (verbose > 0) {
    cout << "  testing   (" << true_expn << ", " << true_digits << ")" << " with n = " << n << endl;
    if (!p) cout << "  got       (" << expn << ", " << buf << ")" << endl;
  }
  pass &= p;
}

bool test_to_digits() {
  cout << "Testing to_digits..." << endl;

  mp_real a;
  int n;
  bool pass = true;

  for (n = 1; n < 15; n++)
    process_to_digits(mp_real(0.0), n, "0", 0, 1, pass);

  for (n = 1; n < 15; n++)
    process_to_digits(mp_real(5.0), n, "5", 0, 1, pass);

  a = 12345.0;
  process_to_digits(a, 1, "1", 4, 1, pass);
  process_to_digits(a, 2, "12", 4, 2, pass);
  process_to_digits(a, 3, "123", 4, 3, pass);
  process_to_digits(a, 4, "1235", 4, 4, pass);
  for (n = 5; n < 15; n++)
    process_to_digits(a, n, "12345", 4, 5, pass);

  a = 1.234999e-40;
  process_to_digits(a, 1, "1", -40, 1, pass);
  process_to_digits(a, 2, "12", -40, 2, pass);
  process_to_digits(a, 3, "123", -40, 3, pass);
  for (n = 4; n <= 6; n++)
    process_to_digits(a, n, "1235", -40, 4, pass);
  for (n = 7; n <= 14; n++)
    process_to_digits(a, n, "1234999", -40, 7, pass);
  process_to_digits(a, 24, "123499900000000000053266", -40, 24, pass);

  a = 9.99999999;
  for (n = 1; n <= 8; n++)
    process_to_digits(a, n, "1", 1, 1, pass);
  for (n = 9; n <= 14; n++)
    process_to_digits(a, n, "999999999", 0, 9, pass);

  return pass;
}

void process_test_read(string in_str, string true_str, bool &pass) {
  mp_real x;
  string str;
  bool p;

  if (x.read(in_str)) {
    str = x.to_string();
    p = (str == true_str);
    if (verbose > 0) {
      cout << "       in = " << in_str << endl;
      cout << "      out = " << str << endl;
      if (!p) cout << "should be = " << true_str << endl;
    }
  } else {
    p = false;
    if (verbose) cout << "read failed." << endl;
  }

  pass &= p;
}

bool test_read() {
  bool pass = true;
  cout << "Testing read..." << endl;

  process_test_read("10 ^ 43 x 1234.5678", "10 ^ 46 x 1.2345678", pass);
  process_test_read("10 ^-43 x 1234.5678", "10 ^ -40 x 1.2345678", pass);
  process_test_read("10 ^43 x-1234.5678", "10 ^ 46 x -1.2345678", pass);
  process_test_read("10 ^-43 x-1234.5678", "10 ^ -40 x -1.2345678", pass);
  process_test_read("10 ^-43 x-1234.56789012345678901234567890", 
                    "10 ^ -40 x -1.2345678901234567890123456789", pass);

  process_test_read("10 ^ 43 x .5678", "10 ^ 42 x 5.678", pass);
  process_test_read("10 ^ 43 X 5678", "10 ^ 46 x 5.678", pass);
  process_test_read("10 ^ 43 x 5678.", "10 ^ 46 x 5.678", pass);

  process_test_read("10", "10 ^ 1 x 1", pass);
  process_test_read("-10", "10 ^ 1 x -1", pass);
  process_test_read("0", "10 ^ 0 x 0", pass);
  process_test_read("-0", "10 ^ 0 x 0", pass);
  process_test_read("12340", "10 ^ 4 x 1.234", pass);
  process_test_read("-12340", "10 ^ 4 x -1.234", pass);
  process_test_read("0.0001234", "10 ^ -4 x 1.234", pass);
  process_test_read("-0.0001234", "10 ^ -4 x -1.234", pass);

  process_test_read("0e0", "10 ^ 0 x 0", pass);
  process_test_read("0e30", "10 ^ 0 x 0", pass);
  process_test_read("0e-30", "10 ^ 0 x 0", pass);
  process_test_read("12.34E456", "10 ^ 457 x 1.234", pass);
  process_test_read("12.34E-456", "10 ^ -455 x 1.234", pass);
  process_test_read("12.34 e 1024", "10 ^ 1025 x 1.234", pass);

  return pass;
}

void process_test(bool pass) {
  if (pass)
    cout << "Passed test." << endl;
  else
    cout << "FAILED test." << endl;
  all_pass &= pass;
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

int main(int argc, char **argv) {
  /* Parse the arguments. */
  char *arg;
  int nr_digits = 100;
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
      verbose++;
    } else {
      cerr << "Unknown flag `" << arg << "'." << endl;
    }
  }

  mp::mp_init(nr_digits + 5);

  process_test(test_binary_io());
  process_test(test_to_digits());
  process_test(test_read());
  
  return all_pass ? 0 : 1;
}
