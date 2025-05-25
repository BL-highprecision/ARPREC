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
#include <cmath> 
#include <arprec/mp_real.h>
#include "small_inline.h"

using std::cerr;
using std::endl;


/*
   Source: Richard Crandall's "Topics in Advanced Scientific Computation",
   p. 82-84.
 
	Alex Kaiser, LBNL, 2/22/2010
 */


mp_real_temp gamma(const mp_real& t) {
	
	// declarations 
	int k; 
	int a, N, sign;
    int working_prec; 
	int ndp; // num digits precision requested. 
	
	ndp = mp::mpgetprec(); 

    working_prec = ceil(1.5 * ndp); 
    // initialize the library to 1.5 times the current precision if necessary
    if(mp::n_digits < working_prec)
        mp::mp_init(working_prec + 10); 
    
    // set the library to 1.5 times the current precision
	mp::mpsetprec(working_prec) ; 
	
    
    a = ceil( 1.25 * ndp / log10( 2.0 * acos(-1.0) ) ) ; 
	
	mp_real arg, sum, runningExp, runningFactorial, temp, x1, x2 ; 
	mp_real denom; 
	// multi-precision constants 
	mp_real rootTwoPi, oneOverRootTwoPi, e, oneOverE ; 
	
	
	// Handle improper arguments.
	if (abs(t) > 1.0e8) {
		cerr << "gamma: agrument is too large" << endl;
		exit(-1); 
	} else if (t == anint(t) && t <= 0.0) {
		cerr << "gamma: invalid negative argument" << endl;
		exit(-1); 
	}
	
	// for testing first handle args greater than 1/2
	// expand with branch later. 
	
	if( t < 0.5 ){
		arg = 1.0 - t; 
		
		// divide by zero trap for later compuation of cosecant
		if( sin(mp_real::_pi * t) == mp_real(0.0) ){
			cerr << "value of argument is too close to a negative integer or zero.\n" << 
					"sin(pi * t) is zero indicating singularity. Increase precision to fix. " << endl; 
			exit(-1); 
		}
	}
	else {
		arg = t; 
		
		// quick exit with factorial if integer
		if(t == aint(t)){
			temp = 1.0; 
			for (k=2; k<t; k++) {
				temp *= k; 
			}
			return temp.toTempAndDestroy();
		}
		
	}
		
	N = a - 1; 
	sign = -1; 
	
	rootTwoPi = sqrt( 2.0 * mp_real::_pi ); 
	oneOverRootTwoPi = mp_real(1.0) / rootTwoPi; 
	
	e = exp(mp_real(1.0) ); 
	oneOverE = mp_real(1.0) / e; 
	runningExp = exp( mp_real(a) ) ; 
	runningFactorial = 1.0; 
	
	sum = 1.0; 
	
	// get summation term
	for(k=1; k <= N; k++){
		sign = -sign; 
		
		// keep (k-1)! term for computing coefficient
		if (k == 1) 
			runningFactorial = 1; 
		else 
			runningFactorial *= k-1; 
		
		runningExp *= oneOverE; // e ^ (a-k). divide by factor of e each iteration 
		
		x1 = mp_real( a - k ); 
		x2 = mp_real(k - 0.5); 
		
		sum += oneOverRootTwoPi * sign * runningExp * pow(x1, x2) / ( runningFactorial * (arg + mp_real(k - 1.0) ))  ; 
	}

    // restore the original precision 
	mp::mpsetprec( ndp ) ;
	
	// compute using the identity
	if( t < 0.5 ){ 
		temp = rootTwoPi * pow(arg + a - 1.0, arg - 0.5) * exp(-arg - a + 1.0) * sum; 
		return mp_real::_pi / (sin(mp_real::_pi * t) * temp) ; 
	}
		
	return rootTwoPi * pow(arg + a - 1.0, arg - 0.5) * exp(-arg - a + 1.0) * sum ;  

}

