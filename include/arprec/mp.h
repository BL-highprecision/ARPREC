/*
 * include/mp.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2001
 *
 */
#ifndef ARPREC_MP_H
#define ARPREC_MP_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
// Above are needed for std::abs, std::min, and std::max.

#include <arprec/arprec_config.h>

#if !(ARPREC_INLINE)
#define inline
#endif

/*
 * This class contains shared variables, macros, and functions that
 * are used by other mp-classes that inherit it.
 */

/* Constants */

/* This gives the index of the first mantissa word in a double MP array.
   The double MP array A is structured as follows:
       Index                   Content
       -----                   -------
       0                    :  the array size
       1                    :  number of mantissa words, with sign
       2                    :  exponent
       3 .. |A[1]|+2        :  mantissa words
       |A[1]|+3 .. |A[1]|+4 :  2 scratch words for rounding purpose
   So, total #words A[0] = |A[1]|+5
*/
#define FST_M       3

/* The following parameters are all that need to be changed in normal usage:

  MPIPL   Initial precision level, in digits.
  MPIOU   Initial output precision level, in digits.
  MPIEP   Log_10 of initial MP epsilon level.
*/
struct ARPREC_API mp_real;

struct ARPREC_API mp {
  // mp_init initializes the ARPREC library with the given precision level.
  // If compute_consts is set to false, it will not compute the constants
  // _pi, _log2, and _log10.  If compute_consts = true and filename is not 
  // NULL, it will attempt to read these constants from the file given by 
  // filename.
  static void mp_init(int new_digits, const char *filename = NULL, bool compute_consts = true);

  static void mp_finalize();

  static void mp_read_constants(const char *filename);
  static void mp_write_constants(const char *filename);

  static void mpsetoutputprec(int num_digits);
  static int mpgetoutputprec()  {return n_output_digits;}
  static void mpsetprec(int num_digits);
  static int mpgetprec();
  static void mpsetprecwords(int num_words);
  static int mpgetprecwords() {return prec_words;}

  static int n_digits;
  static int n_mantissa;
  static int n_words;
  static int n_output_digits;

  static int fmpwds5; // the static word size used in Fortran 90 wrapper.

  static int debug_level; // Defaults to 0.
  static int debug_words; // Words in debug output.  Defaults to 22.
  static int error_no;    // Error indicator.  Initially 99.

  // Note: these are approximate in the sense that the result is computed
  //       to extra word, and then rounded to the given direction.
  enum rounding_mode {
    round_to_zero = 0, 
    round_to_nearest = 1, 
    round_to_infinity = 2
  };

  // ARPREC rounding mode. Defaults to round_to_zero.
  static enum rounding_mode round_dir;
	
  // Array of error options. Initial value: 2.
  static int MPKER[79+1];
  
  // array pointers for mpmulx, other high precision routines.
  // They store roots of unity (complex).
  static double *mpuu1;
  static double *mpuu2;
  static double *mpuu3;

  // used to store old FPU control word, if any.
  // set and used by mp_init and mp_finalize.
  static unsigned int old_cw;

  // Current precision level, in words. Initial value: 16.
  static int prec_words; // n_mantissa + 1

  // Depends on the system word size.
  //  On IEEE systems and most other 32 bit systems, set to 4096.D0
  static const double mpbbx;
  static const double mpbdx; 
  static const double mprdx; 
  static const double mpbx2; 
  static const double mprbx; 
  static const double mprx2; 
  static const double mprxx;
  
  // Depends on the system word size.
  // On IEEE systems and most other 64 bit systems, set to 52.
  static const int mpnbt;
	
  // Depends on the system word size.
  // On IEEE systems and most other 64 bit systems, set to 64.
  static const int mpnpr;
	
  // Depends on the system word size.
  // On IEEE systems and most other 32 bit systems, set to 7.
  // On Cray system, set to 8.
  static const int mpmcrx;
	
  // This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  // Value of 32 appears to work well on most systems.
  static const int mpnrow;
	
  // This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  // Value of 3 appears to work well on most systems. 
  static const int mpnsp1; 
	
  // This is spacing parameters to avoid bank and cache conflicts in the FFT routines.
  // Value of 17 appears to work well on most systems.
  static const int mpnsp2;

  static const double digits_per_word;

  //  This routine terminates execution. 
  static void mpabrt();

  static const double _d_nan;
};

#endif
/* ARPREC_MP_H */
