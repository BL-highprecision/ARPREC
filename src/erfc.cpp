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
 *  Handles the error: MPID
 *  if it is greater=4
 */
static void erfMPIDErrorHandle(const mp_real& arg) {
  if(mp::debug_level >= 4) {
    int no = static_cast<int>(arg.mpr[0]);
    no = std::min(no,mp::debug_words) + 2;
    cerr << "arg.mpr[0] = " << arg.mpr[0] << " no = " << no << endl;
  }
}

/*
 * calculates erfc(t) = 1 - erf(t)
 *
 */
mp_real  erfc(const mp_real& arg)
{
    //Static Variables that are used to make the function
    //run faster

    //Initialization of Static Variables

    //current size     Previous Prec
    static int size=0, pPrec = -1;
    static double eps[8] = {8,1,4-mp::n_words,1,0,0,0,0};
    static mp_real eps_m = mp_real(eps);
    //Table containing exp(-k^2 * alpha^2)
    //but k=1 one starts at index 0
    double alpha = DBL_MAX, alpha2 = 0.0;
    static const double dpi = dble(mp_real::_pi),
	dlog10 = dble(mp_real::_log10),
	dlog2  = dble(mp_real::_log2);
    static mp_real *epow;

    //Error Checking
    erfMPIDErrorHandle(arg);
    if(mp::error_no!=0){
	cerr << "error_no Error " << mp::error_no << endl;
	return mp_real(0.0);
    }

    //Exit Conditions/ Result Checks
    if(arg <= 0.0)
	return (arg == 0.0) ? mp_real(1.0) : 2.0 - erfc(-arg);
    //Arg is too big virtually zero for our purpose
    //This is a hidden bottle neck on precesion
    if(arg > 10000.0)
	return mp_real(0.0);

    mp_real sk0, sk1, sk2, sk3, sk4, sk5, sk6;
    
    int n;
    const int currentPrec=mp::mpgetprec();
    double d1,d2;

    //converting arg to double format
    mp_real::mpmdc(arg, d1, n, mp::prec_words);
    d1 = ldexp(d1, n);

    d2 = dpi/d1;

    /*
     *  ePow table is going to be recalculated
     *  On first call or increase in working precession
     *  or a need to increse working precession
     */
    if( d2 < alpha || pPrec < currentPrec){
	pPrec = currentPrec;  //setting the Precesion

	//deleting all the allocated memory
	if(epow){
            delete [] epow;
            epow = NULL;
	}


	//Calculating the table  of exp(-k^2*alpha^2)

	//inits for the talble calculation
	double sqrtplog10 = sqrt(currentPrec*dlog10);
	d1 = dpi/sqrtplog10;
	if(d1 > d2) d1 = d2;
	// Mulltiply d1 (new alpha) by 0.95 (so we won't need to recalculate
	// so often), then round to some nice 6-bit rational.
	d1 *= .95;

	n = abs(static_cast<int>(log(d1)/dlog2))+1;
        alpha = ldexp(anint(ldexp(d1, n+6)), -n-6);
	size = static_cast<int>(sqrtplog10/alpha) + 1;

	//Make sure that (alpha * ntab)^2 can be represented exactly in DP.
	//I don't think this will ever be a problem, but check just in case.
	d2 = 2*(6.0+log(static_cast<double>(size))/dlog2);
	if(d2 > 53.0){
	    cerr <<"erfc: error; contact author" << endl;
	    exit(0);
	}

	//Allocation of needed memory
        epow = new mp_real[size];
	if(!epow){
	    cerr << "QuadErf::erfc: Not enough Memory for e^() table" << endl;
	    exit(-1);
	}

	//t1 = - alpha ** 2
	alpha2 = alpha*alpha;
        mp_real::mpdmc(-alpha2, 0, sk0, mp::prec_words);

	//t2 = exp (t1)
        mp_real::mpexpx(sk0, mp_real::_pi, mp_real::_log2, sk1);

	// t3 = t2 ** 2
        mp_real::mpsqx(sk1, sk2, mp::prec_words);

	//t4 = 1.d0
        sk3 = 1.0;

	for(int i=0; i< size; i++){
	    //t4 = t2 * t4
            mp_real::mpmulx(sk1, sk3, sk4, mp::prec_words);
            mp_real::mpeq(sk4, sk3, mp::prec_words);
	    
	    //ePow[i] = t4
            mp_real::mpeq(sk3, epow[i], mp::prec_words);

	    //t2 = t2 * t3  (oddSequence)
            mp_real::mpmulx(sk1, sk2, sk4, mp::prec_words);
            mp_real::mpeq(sk4, sk1, mp::prec_words);

	}
    }

    //t1 = 0.d0
    sk0 = 0.0;

    //t2 = t ** 2
    mp_real::mpsqx(arg, sk1, mp::prec_words);

    //t3 = exp (-t2)
    sk4 = -sk1;
    mp_real::mpexpx(sk4, mp_real::_pi, mp_real::_log2, sk2);
    
    for(int k=0; k<size; k++){
	//t5 = ePow[k] / (k2 ** 2 * alpha ** 2 + t2)
	int k2 = k+1;
        mp_real::mpdmc((k2*k2)*alpha2, 0, sk4, mp::prec_words);
        mp_real::mpadd(sk4, sk1, sk5, mp::prec_words);
        mp_real::mpdivx(epow[k], sk5, sk4, mp::prec_words);

	//t1 = t1 + t5
        mp_real::mpadd(sk0, sk4, sk5, mp::prec_words);
        mp_real::mpeq(sk5, sk0, mp::prec_words);

	//if (abs (t5) < t4) break;
        sk5 = abs(sk4);
        if (sk5 < eps_m) break;
    }
    

    //erfc = t3 * alpha * t / PI * (1.d0 / t2 + 2.d0 * t1) &
    //       + 2.d0 / (1.d0 - exp (2.d0 * PI * t / alpha))
    mp_real::mpmuld(sk2, alpha, 0, sk3, mp::prec_words);
    mp_real::mpmulx(sk3, arg, sk4, mp::prec_words);
    mp_real::mpdivx(sk4, mp_real::_pi, sk3, mp::prec_words);
    mp_real::mpdivx(mp_real(1.0), sk1, sk4, mp::prec_words);
    mp_real::mpmuld(sk0, 2.0, 0, sk5, mp::prec_words);
    mp_real::mpadd(sk4, sk5, sk6, mp::prec_words);
    mp_real::mpmulx(sk3, sk6, sk4, mp::prec_words);
    mp_real::mpmuld(mp_real::_pi, 2.0, 0, sk3, mp::prec_words);
    mp_real::mpmulx(sk3, arg, sk5, mp::prec_words);
    mp_real::mpdivd(sk5, alpha, 0, sk6, mp::prec_words);
    mp_real::mpexpx(sk6, mp_real::_pi, mp_real::_log2, sk3);
    mp_real::mpsub(mp_real(1.0), sk3, sk5, mp::prec_words);
    mp_real::mpdmc(2.0, 0, sk6, mp::prec_words);
    mp_real::mpdiv(sk6, sk5, sk3, mp::prec_words);

    mp_real ans;
    mp_real::mpadd(sk4, sk3, ans, mp::prec_words);
    erfMPIDErrorHandle(arg);
    return ans;
}

