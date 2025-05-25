/*
 * src/binary_io.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2006
 *
 * Handles binary I/O of mp_real
 */
#include <cstdlib>
#include <limits>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

// Dumps the internal structure of mp_real number A to the stream 
// out.  Beware that the output formatting is machine dependent due
// to endianness.
void mp_real::write_binary(std::ostream& out) {
  int len;

  // write the significant part
  len = static_cast<int>(std::abs(mpr[1])) + FST_M;
  out.write(reinterpret_cast<char *>(mpr), sizeof(double) * len);

  // write out the zeros
  len = static_cast<int>(mpr[0]) - len;
  double zero = 0.0;
  for (int i = 0; i < len; i++)
    out.write(reinterpret_cast<char *>(&zero), sizeof(double));
}

// Reads an mp_real number from the stream in.  The input should
// be a binary array of doubles.  If the input is larger than the
// space allocated for *this, it is truncated.  Beware the formatting
// is machine dependent due to endianness.
void mp_real::read_binary(std::istream &in) {
  double v;
  int len;

  // read in the array length
  in.read(reinterpret_cast<char *> (&v), sizeof(double));
  if (v < FST_M + 2) {
    cerr << "*** mp_real::read_binary: array length too short." << endl;
    mpabrt();
  } else if (v > std::numeric_limits<int>::max()) {
    cerr << "*** mp_real::read_binary: array length too long." << endl;
    mpabrt();
  }

  // read in min ( space in mpr, array length in file )
  len = static_cast<int>(std::min(v, mpr[0]));
  in.read(reinterpret_cast<char *>(&mpr[1]), sizeof(double)*(len - 1));
  if (in.eof()) {
    cerr << "*** mp_real::read_binary: unexpected EOF." << endl;
    mpabrt();
  }

  // truncate the extraneous entries
  len = static_cast<int>(v) - len;
  in.ignore(len * sizeof(double));
  if (in.eof()) {
    cerr << "*** mp_real::read_binary: unexpected EOF." << endl;
    mpabrt();
  }
  
  mpr[1] = sign(1.0, mpr[1]) * std::min(mpr[1], mpr[0]-FST_M-2);
}

