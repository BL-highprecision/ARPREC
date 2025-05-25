/*
 * src/mpreal_friends.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 */
#include <cfloat>
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;

/*
 * The Error Function
 * erf = Int_0^a 2/Sqrt(pi) * e^(-arg^2)
 */
mp_real_temp erf(const mp_real& arg)
{
    if ( abs(arg) > 1.0e-4 ) {
	return (1.0 - erfc(arg));
    } else {
	double ds;
	mp_real eps, t0, t1, t2, t3, t4;
	int nwords = mp::prec_words, i;

	eps = pow(mp_real(0.5), mp::mpnbt * (nwords + 1));
	t0 = t1 = arg;
	t2 = pow(arg, 2);
	t3 = ds = 1.0;

	for (i = 1; i < 1000000000; ++i) {
	    ds = -ds;
	    t3 = t3 * double(i);
	    t1 = t1 * t2;
	    t4 = ds * t1 / (double(2 * i + 1) * t3);
	    t0 = t0 + t4;
	    if ( abs(t4) < eps ) break;
	}
	if ( i == 1000000000 ) {
	    cerr << "erf: loop end error" << endl;
	    t0 = 0.0;
	}
	return (2.0 / sqrt (mp_real::_pi) * t0);
    }
}

