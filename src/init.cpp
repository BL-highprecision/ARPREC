/*
 * src/mpreal2.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2002
 *
 */
#include <fstream>
#include <limits>
#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>
#include <arprec/mp_int.h>
#include <arprec/fpu.h>

#if !(ARPREC_INLINE)
#include <arprec/mp_inline.h>
#include <arprec/mp_complex_inline.h>
#include <arprec/mp_int_inline.h>
#endif

using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;

// initializing static variables
const double mp::mpbbx = 16777216.0;//33554432.0;  // 2^25
const double mp::mpbdx = mp::mpbbx * mp::mpbbx;// 2^50
const double mp::mpbx2 = mp::mpbdx * mp::mpbdx; 
const double mp::mprbx = 1.0 / mp::mpbbx;
const double mp::mprdx = 1.0 / mp::mpbdx; // 2^(-50)
const double mp::mprx2 = mp::mprdx * mp::mprdx; 
const double mp::mprxx = 16.0 * mp::mprx2;

const double mp::_d_nan = std::numeric_limits<double>::quiet_NaN();

const double mp::digits_per_word = mp::mpnbt * 3.01029995663981195214e-01; // log(2) / log(10);

const int mp::mpnbt = 48, //bits ber word
  mp::mpnpr = 16,//4, //half the number of words added before forcing carry
  
  mp::mpmcrx = 7, // advanced routines start at roughly 2^mpmcrx 
  mp::mpnrow = 16,//32
  mp::mpnsp1 = 3, //3
  mp::mpnsp2 = 9; //17


enum mp::rounding_mode mp::round_dir = mp::round_to_nearest;
int mp::debug_level=0,
  mp::debug_words=22,
  mp::error_no=99,
  mp::MPKER[79+1] = {0,2,2,2,2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2,2,2,2,2,
                            2,2,2,2,2,2};

int mp::n_mantissa;
int mp::prec_words,    // all these get set later in mp_init
  mp::n_output_digits,
  mp::n_words,
  mp::n_digits;

int mp::fmpwds5; // the static word size used in Fortran 90 wrapper,
                 // fmpwds5 >= n_words, and is set by c_mpinit() in c_mp.cc.

mp_real mp_real::_pi(0.0, 0);
mp_real mp_real::_log2(0.0, 0);
mp_real mp_real::_log10(0.0, 0);
mp_real mp_real::_eps(0.0, 0);
double *mp::mpuu1 = 0;
double *mp::mpuu2 = 0;
double *mp::mpuu3 = 0;

unsigned int mp::old_cw = 0;

static void init_constants();
static void init_eps();

void mp::mp_init(int new_digits, const char *filename, bool compute_consts) {
  fpu_fix_start(&mp::old_cw);

  const double digits_per_word = mp::mpnbt * log(2.0)/log(10.0);

  n_digits = new_digits;
  n_mantissa = int((n_digits-1) / digits_per_word + 2);
  prec_words = n_mantissa;
  n_words = prec_words + 5;

  mp::n_output_digits = new_digits;

  //This lets other procedures know that the library is ready.
  error_no = 0;
  if(prec_words > (1<<(mpmcrx-1))) 
    mp_real::mpinix(prec_words+8);

  mp_real::_eps.reallocate(7);
  init_eps();

  if (compute_consts) {
    if (filename)
      mp_read_constants(filename);
    else
      init_constants();
  }
}

void mp::mp_finalize() {
  fpu_fix_end(&mp::old_cw);
  if (mp::mpuu1) {
    delete [] mp::mpuu1;
    mp::mpuu1 = NULL;
  }
  if (mp::mpuu2) {
    delete [] mp::mpuu2;
    mp::mpuu2 = NULL;
  }
  if (mp::mpuu3) {
    delete [] mp::mpuu3;
    mp::mpuu3 = NULL;
  }
}

void mp::mpsetoutputprec(int num_digits) {
  n_output_digits = std::min(n_digits-2, std::max(1, num_digits));
  if(n_output_digits != num_digits) {
    cerr << "Request for output of " << num_digits 
         << " did not succeed." << endl;  
    cerr << "MPINIT must first be called with at least "
         << num_digits+2 << " digits of precision." << endl;
    cerr << "Defaulting to output of " << n_output_digits << " digits." << endl;
  }
}

void mp::mpsetprec(int num_digits) {
  mpsetprecwords(static_cast<int>( (num_digits-1) / digits_per_word ) + 2);
}

int mp::mpgetprec() {
  return static_cast<int>( (prec_words-1) * digits_per_word) + 1; 
}

void mp::mpsetprecwords(int num_words) {
  prec_words = std::max(0, std::min(n_mantissa+1, num_words));
  n_words = prec_words + 5;
}

void mp::mpabrt()
{
  if (error_no == 99) 
    cerr << "*** The ARPREC library has not been initialized." << endl;
  else
    cerr << "*** mpabrt: execution terminated, error code =" << error_no << endl;
  exit(-1);
}

static void init_pi() {
  mp::prec_words += 4;
  mp_real::_pi.reallocate(mp::prec_words + 6);
  mp_real::mppix(mp_real::_pi);
  mp::prec_words -= 4;
}

static void init_log2() {
  mp_real::prec_words+=4;
  int prec_words = mp_real::prec_words;
  mp_real::_log2.reallocate(prec_words + 5);
  mp_real::_log2 = 3.0;
  mp_real t2(2.0, 6);
  // fourth optional "hidden" nit argument set to zero for greater accuracy
  // users should not call mplog with four arguments, only 3.
  if(prec_words < (1<<(mp::mpmcrx-2)))
    mp_real::mplog(t2, mp_real::_log2, mp_real::_log2, prec_words, 0);
  else
    mp_real::mplogx(t2, mp_real::_pi, mp_real::_log2, mp_real::_log2, prec_words);
  mp_real::prec_words-=4;
}

static void init_log10() {
  mp_real::prec_words+=4;
  int prec_words = mp_real::prec_words;
  mp_real::_log10.reallocate(prec_words + 5);
  mp_real t2(10.0, 8);
  mp_real::mplogx(t2, mp_real::_pi, mp_real::_log2, mp_real::_log10, prec_words);
  mp_real::prec_words-=4;
}

static void init_constants() {
  init_pi();
  init_log2();
  init_log10();
}

static void init_eps() {
  mp_real t2(10.0, 6);
  mp_real::prec_words+=4;
  int prec_words = mp::prec_words;
  mp_real::mpnpwx(t2, -mp::n_digits, mp_real::_eps, prec_words);
  mp_real::prec_words-=4;
}

void mp::mp_read_constants(const char *filename) {
  ifstream infile(filename);
  if(!infile) {
    cerr << "Could not open MP initialization file "<< filename << endl;
    mpabrt();
  }
  int temp;
  infile >> temp;
  if(temp < n_digits) {
    cerr << "MP Initialization file incorrect or does not "
         << "have sufficient precision." << endl;
    mpabrt();
  }
  infile.ignore(); //get rid of newline.
  mp_real::_pi.read_binary(infile);
  mp_real::_log2.read_binary(infile);
  mp_real::_log10.read_binary(infile);

  // sanity check read constants
  double t1;
  int n1;
  mp_real::mpmdc(mp_real::_pi, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - 3.141592653589793) > mprx2) {
    cerr << "*** MPINIT: Pi is wrong in file " << filename << endl;
    mpabrt();
  }
  mp_real::mpmdc(mp_real::_log2, t1, n1, prec_words);
  if(n1 != -mpnbt || (std::abs(t1 * mprdx - 0.693147180559945309) > mprx2)) {
    cerr << "*** MPINIT: Log(2) is wrong in file " << filename << endl;
    mpabrt();
  }
  mp_real::mpmdc(mp_real::_log10, t1, n1, prec_words);
  if(n1 != 0 || std::abs(t1 - 2.3025850929940459) > mprx2) {
    cerr << "*** MPINIT: Log(10) is wrong in file " << filename << endl;
    mpabrt();
  }
}

void mp::mp_write_constants(const char *filename) {
  ofstream outfile(filename);
  if(!outfile) {
    cerr << "Cannot output initialization file " << filename << endl;
    return;
  }
  outfile << mp::n_digits << endl;
  mp_real::_pi.write_binary(outfile);
  mp_real::_log2.write_binary(outfile);
  mp_real::_log10.write_binary(outfile);
  outfile.close();
}
