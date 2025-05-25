/*
 * src/debug.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2007
 *
 */
#include <iostream>
#include <cmath>
#include <arprec/mp_real.h>

using std::cerr;
using std::endl;

void mp_real::print_short(const char *name) const {
  int nw = static_cast<int>(mpr[0]);
  int nm = static_cast<int>(std::abs(mpr[1]));
  int e  = static_cast<int>(mpr[2]);

  cerr << name << " = " << "[" << nw << ", " << nm << "] ";
  cerr << "[ " << e << "]";
  if (nm <= 0)
    cerr << "0";
  else
    cerr << mpr[3];
}

void mp_real::print_mpreal(const char* str, 
    const mp_real& mpr, std::ostream &os) {
  int nm, nw, w;
  
  nw = static_cast<int>(mpr[0]);
  nm = std::abs(static_cast<int>(mpr[1]));
  os << str << "  nw = " << nw << endl;
  os << str << "  nm = " << nm << endl;
  os << str << " exp = " << mpr[2] << endl;
  os.precision(20);
  for(w = FST_M; w < std::abs(mpr[1]) + FST_M; ++w)
    os << "  w = " << w << "  " << mpr[w] << endl;
  os.precision(0);
}

