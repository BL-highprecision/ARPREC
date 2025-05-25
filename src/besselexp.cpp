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
 * This Evalutes the function BesselI (0, t).
 * July 9th, 2004
 */
mp_real_temp besselexp(const mp_real& t)
{
    const mp_real c10(10.0);
    int currentPrec = mp::mpgetprec();
    const mp_real eps(pow(c10, -currentPrec - 10));
    mp_real sum(1.0), t1(1.0);

    //Select eithter the direct or the asymptotic series
    if(.85*t < currentPrec)
    {
	mp_real t2 = sqr(t);
	for(int i=1; i<=10000000; i++)
	{
	    t1 *= 0.25 * t2 / sqr(mp_real(i));
	    if(t1 < eps)
		return sum/exp(t);
	    sum += t1;
	}
	cerr << "besselexp: loop overflow 1" << endl;
    }else
    {
	for(int i=1; i<=10000000; i++)
	{
	    const mp_real t2(t1);
	    t1 *= sqr(mp_real(2*i-1.0))/((i*8.0)*t);
	    sum += t1;

	    if(t1 < eps)
		return sum/sqrt(2.0*mp_real::_pi*t);

	    if(t1 > t2){
		cerr << "besselexp: t1 > t2; t = " << t << endl;
		return mp_real(0.0).toTempAndDestroy();
	    }
	}
	cerr << "besselexp: loop overflow 2" << endl;
    }
    return mp_real(0.0).toTempAndDestroy();
}

