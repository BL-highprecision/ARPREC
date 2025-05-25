/*
 * quad-erf.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 7, 2004
 * Modified:      June 27, 2004
 * Version:       1.0
 */

//Implementation of QuadErf
#include <cmath>
#include <cfloat>
#include <iostream>
#include "quad-erf.h"
#include "util.h"

using std::cout;
using std::endl;

//Static Variable initialization
int      QuadErf::_phasesMax    = 0;
long int QuadErf::_AbWtSize     = 0;
long int QuadErf::_AbWtMax      = 0;
long int QuadErf::_quadErfCount  = 0;
mp_real* QuadErf::_abscissas    = NULL;
mp_real* QuadErf::_weight       = NULL;
int      QuadErf::_precWord1    = 0;
int      QuadErf::_precWord2    = 0;
int      QuadErf::_sneps1       = 0;
int      QuadErf::_sneps2       = 0;


/*
 *  ABSTRACT FUNCTION
 *  Resizes the Abscissas and weight arrays
 *
 *  This subroutine initializes the quadrature arays xk and wk using the
 *  function x(t) = tanh (pi/2*sinh(t)).  The argument nq2 is the space
 *  allocated for wk and xk in the calling program.  By default it is set to 
 *  12 * 2^nq1.  Increase nq2 if directed by a message produced below.
 *  Upon completion, wk(-1) = nq1, and xk(-1) = n, the maximum space parameter
 *  for these arrays.  In other words, the arrays occupy (wk(i), i = -1 to n)
 *  and (xk(i), i = -1 to n), where n = xk(-1).   The array x_k contain 1 minus
 *  the abscissas; the wk array are the weights at these abscissas.
 */
void QuadErf::reSizeAbWt(long int AbWtSize)
{
    if(_phasesMax < _phases)
	_phasesMax = _phases;

    if(AbWtSize < 12*2*2)
    {
	AbWtSize = 12;
	//AbWtSize = 8;
	for(int i=0; i<_phasesMax; i++)
	    AbWtSize *= 2;
    }
    //Need to recompute everything as precesion has been increased
    if(_precWord1 < _precWd1){
	_AbWtSize = _AbWtMax=  0;
	_precWord1 = _precWd1;
	_sneps1    = _neps1;
	if(_abscissas)
	{
	    delete [] _abscissas;
	    delete [] _weight;
	    _abscissas = _weight = NULL;
	}
    }

    //There is no need to inc to the given Size
    if((_precWord2 >= _precWd2  || _AbWtSize == _AbWtMax) && 
       _AbWtSize >= AbWtSize)
	return;

    if(_precWord2 < _precWd2){
	_precWord2 = _precWd2;
	_sneps2    = _neps2;
    }

    mp::mpsetprecwords(_precWord1);

    //increasing the table Maxsize only if it is bigger than the 
    //current table Maxsize
    if(_AbWtMax < AbWtSize)
    {
	_AbWtMax = AbWtSize;
	mp_real* ta = new mp_real[_AbWtMax]; //abscissas
	mp_real* tw = new mp_real[_AbWtMax]; //weight

	if(!(ta && tw))
	{
	    cout << "QuadErf::reSizeAbWt: Not enough memory" <<endl;
	    exit(0);
	}

	//Copying the values from the old array to the new
	for(long int i = 0; i < _AbWtSize; i++)
	{
	    ta[i] = _abscissas[i];
	    tw[i] = _weight[i];
	}

	//Making sure that arrays are not NULL
	//if one is null other must be null as they are the same size
	if(_abscissas)
	{
	    delete [] _abscissas;
	    delete [] _weight;
	}
	_abscissas = ta;
	_weight = tw;
    }

    if (_ndebug >= 1)
	cout <<"QuadErf::reSizeAbWt: Error Function quadrature initialization"
	     <<endl;
    
    /*
     *  Some inits before actual computation of tables begin
     *
     */
    mp::mpsetprecwords(_precWord2);
    mp_real eps2 = pow(mp_real(10.0), _sneps2);
    mp::mpsetprecwords(_precWord1);
    mp_real p2   = .5*mp_real::_pi;
    mp_real spi  = sqrt(4.0/mp_real::_pi);
    mp_real h    = pow(mp_real(.5), _phasesMax-2);
#define IPRINT 1000
    //computation for the table
    for(long int i = _AbWtSize; i < _AbWtMax; i++)
    {
	mp_real t1(h*(double)i);
	_abscissas[i] = erfc(t1);
	_weight[i] = spi/exp(sqr(t1));
	
	if(_ndebug >=2 && i%IPRINT == 0){
	    cout << "\t\t" << i << "\t" << _AbWtMax << endl;
	}

	if (isThereErrorMPIER(_ierror, "incAbWtArray: Table space parameter "
			      "is too small"))
	    return;
	
	if(_weight[i] < eps2)
	{
	    _AbWtSize = i;
	    if(_ndebug >= 2){
		cout<<"incAbWtArray: Table spaced used ="<<_AbWtSize<<endl;
		cout << "\t\t" << i << "\t" << _AbWtMax << endl;
	    }
	    return;
	}
    }
#undef IPRINT
    _AbWtSize = _AbWtMax;

    cout << "incAbWtArray: Table space parameter is too small; value ="
	 << _AbWtSize << endl;
    _nerror = 91;
  }



/*
 *  This routine computes the integral of the function in fun on the interval
 *  [x1, x2], with up to nq1 iterations, with a target tolerance of 10^epsilon1.
 *  xk and wk are precomputed tables of 1 - abscissas and weights.  The function
 *  fun is not evaluated at x = x1 or x2.  The array _abscissas contain 1 minus
 *  the abscissas; the _weight array are the weights at these abscissas.
 */
mp_real QuadErf::integrate(mp_real func(const mp_real &x),
			   const mp_real &x1, const mp_real &x2)
{
    static const int izx = 5;
    int outprec = mp::mpgetoutputprec();
    mp::mpsetoutputprec(56);

    //inits
    mp::mpsetprecwords(_precWd2);
    mp_real ax = .5*(x2-x1);
    mp_real bx = .5*(x2+x1);
    mp::mpsetprecwords(_precWd1);
    const mp_real c10(10.0);
    const mp_real ceps1(pow(mp_real(10.0), _neps1));
    const mp_real ceps2(pow(mp_real(10.0), _neps2));
    mp_real sum(0.0),s1(0.0),s2(0.0),s3(0.0),h(4.0),
	fm1(0.0),fm2(0.0);
    mp_real err(0.0), t1(0.0), t2(0.0), tw1(0.0), tw2(0.0);

    int phases = _phasesMax+1, k1, k2, iz1,iz2;
    int *ip = new int[phases];
    //Creating a table to minimize computation
    ip[0] = 1;
    for(int i=1; i < phases; i++)
	ip[i] = 2*ip[i-1];

    for(long int k = 1; k <= _phases; k++)
    {
	h *= .5;
	s3 = s2;
	s2 = s1;
	fm1 = 0.0;
	fm2 = 0.0;
	k1  = ip[_phasesMax-k];
	k2  = ip[_phasesMax-k+1];
	
	iz1 = iz2 = 0;

	//Evaluate function at level k in x, avoiding unnecessary computation.
	for(long int i = 0; i < _AbWtSize; i += k1)
	{
	    if(i%k2 != 0 || k == 1)
	    {
		mp::mpsetprecwords(_precWd2);
		mp_real xt1(1.0 - _abscissas[i]);
		mp_real xx1(bx - ax*xt1);
		mp_real xx2(bx + ax*xt1);
		bool log1 = xx1 > x1 && iz1 < izx;
		bool log2 = xx2 < x2 && iz2 < izx;
		mp::mpsetprecwords(_precWd1);

		if(!log1 && !log2)
		    break;

		if(log1)
		{
		    t1 = func(xx1);
		    if(isThereErrorMPIER(_ierror) || _nerror > 0)
		    {
			if(_ierror > 0) _nerror = 100 + _ierror;
			cout << "QuadErf::integrate: Error in quadrature"
			    "calculation; code =" << _nerror << endl;
			mp::mpsetoutputprec(outprec);
			delete ip;
			return mp_real(0.0);
		    }
		    tw1 = t1*_weight[i];
		    if(abs(tw1) < ceps1) iz1++;
		    else iz1 = 0;
		}else
		{
		    t1  = 0.0;
		    tw1 = 0.0;
		}

		if(i > 0 && log2)
		{
		    t2 = func(xx2);
		    if(isThereErrorMPIER(_ierror) || _nerror > 0)
		    {
			if(_ierror > 0) _nerror = 100 + _ierror;
			cout << "QuadErf::integrate: Error in quadrature"
			    "calculation; code =" << _nerror << endl;
			mp::mpsetoutputprec(outprec);
			delete ip;
			return mp_real(0.0);
		    }
		    tw2 = t2*_weight[i];
		    if(abs(tw2) < ceps1) iz2++;
		    else iz2 = 0;
		}else
		{
		    t2  = 0.0;
		    tw2 = 0.0;
		}

		sum += tw1 + tw2;
		tw1 = abs(tw1);
		tw2 = abs(tw2);
		t1  = abs(t1);
		t2  = abs(t2);

		fm1 = std::max(fm1, tw1);
		fm1 = std::max(fm1, tw2);
		fm2 = std::max(fm2, t1);
		fm2 = std::max(fm2, t2);

	    }
	}


	//Intergral = s1 AND Error estimation calculation follows
	s1 = ax*h*sum;
	mp_real eps1(fm1*ceps1);
	mp_real eps2(fm2*ceps2);
	double d1 = dplog10q(abs(s1 - s2));
	double d2 = dplog10q(abs(s1 - s3));
	double d3 = dplog10q(eps1) - 1;
	double d4 = dplog10q(eps2) - 1;

	if(k <= 2)
	    err = 1.0;
	else if(d1 == -9999.0)
	    err = 0.0;
	else
	{
	    double val1 = d1*d1/d2;
	    double val2 = 2*d1;
	    double max = std::max(val1,val2);
	    max = std::max(max,d3);
	    max = std::max(max,d4);
	    err = pow(c10, nint(std::min(0.0,max)));
	}

	if(_ndebug >= 2)
	    cout << "quadErf: Iteration=" << k << " of " << _phases
		 << "; est error = 10^" << nint(dble(dplog10q(abs(err))))
		 << "; approx value = " << s1 << endl;

	if(k > 3 && err < eps1){
	    mp::mpsetoutputprec(outprec);
	    delete ip;
	    return s1;
	}
	
	if(k >= 3 && err < eps2)
	{
	    cout << "QuadErf::integrate: Estimated error = 10^" 
		 << nint(dble(dplog10q(abs(err)))) << "\n"
		 << "Increase secondary prec (Ndigits2) for greater accuracy. "
		 << "Current value =" << -_neps2 << endl;
	    mp::mpsetoutputprec(outprec);
	    delete ip;
	    return s1;
	}
    }

    cout << "QuadErf::integrate: Estimated error = 10^"
	 << nint(dble(dplog10q(abs(err)))) << "\n"
	 << "Increae QuadLevel for greater accuracy. "
	 << "Current value =" << _phases << endl;

    mp::mpsetoutputprec(outprec);
    delete ip;
    return s1;
}


/*
 *  Handles the error: MPID
 *  if it is greater=4
 *
static bool erfMPIDErrorHandle(const mp_real& arg)
{
    if(mp::debug_level >= 4){
	int no = (int)arg.mpr[0];
	no = std::min(no,mp::debug_words) + 2;
	cout <<"arg.mpr[0]="<< arg.mpr[0]<< " no=" << no << endl;
	return true;
    }
    return false;
}
*/


/*
 * calculates erfc(t) = 1 - erf(t)
 *
 *
mp_real erfc(const mp_real& arg)
{
    //Static Variables that are used to make the function
    //run faster

    //Initialization of Static Variables

    //current size     Previous Prec
    static int size=0, pPrec = -1;
    static double eps[8] = {8,1,4-mp::mp5,1,0,0,0,0};
    static const double one[8] = {8,1,0,1,0,0,0,0};
    //Table containing exp(-k^2 * alpha^2)
    //but k=1 one starts at index 0
    static double **ePow = NULL,
	alpha(DBL_MAX),alpha2(0.0);
    static const double dpi = dble(mp_real::_pi),
	dlog10 = dble(mp::mpl10),
	dlog2  = dble(mp::mpl02);

    //Error Checking
    erfMPIDErrorHandle(arg);
    if(isThereErrorMPIER("error_no Error"))
	return mp_real(0.0);

    //Exit Conditions/ Result Checks
    if(arg <= 0.0)
	return (arg == 0.0) ? mp_real(1.0) : 2.0 - erfc(-arg);
    //Arg is too big virtually zero for our purpose
    //This is a hidden bottle neck on precesion
    if(arg > 10000.0)
	return mp_real(0.0);


    //Initialization of variables
    const double * t = arg.mpr;
    double *sk0 = new double[mp::mp5],
	*sk1 = new double[mp::mp5],
	*sk2 = new double[mp::mp5],
	*sk3 = new double[mp::mp5],
	*sk4 = new double[mp::mp5],
	*sk5 = new double[mp::mp5],
	*sk6 = new double[mp::mp5];
 
    sk0[0]=sk1[0]=sk2[0]=sk3[0]=sk4[0]=sk5[0]=sk6[0]=mp::mp5;
    
    int n;
    const int currentPrec=mp::mpgetprec();
    double d1,d2;

    //converting arg to double format
    c_mpmdc(t,&d1,&n);
    d1 *= pow(2.0,n);

    d2 = dpi/d1;

    // 
    //   ePow table is going to be recalculated
    //   On first call or increase in working precession
    //   or a need to increse working precession
    // 
    if( d2 < alpha || pPrec < currentPrec){
	pPrec = currentPrec;  //setting the Precesion


	//deleting all the allocated memory
	if(ePow){
	    for(int i = 0; i < size; i++)
		delete [] ePow[i];
	    delete [] ePow;
	}


	//Calculating the table  of exp(-k^2*alpha^2)

	//inits for the talble calculation
	double sqrtplog10 = sqrt(currentPrec*dlog10);
	d1 = dpi/sqrtplog10;
	if(d1 > d2) d1 = d2;
	// Mulltiply d1 (new alpha) by 0.95 (so we won't need to recalculate
	// so often), then round to some nice 6-bit rational.
	d1 *= .95;

	n = abs((int)(log(d1)/dlog2))+1;
	alpha = pow(.5,n+6)*anint(d1*pow(2.0,n+6));
	size = (int)(sqrtplog10/alpha) + 1;

	//Make sure that (alpha * ntab)^2 can be represented exactly in DP.
	//I don't think this will ever be a problem, but check just in case.
	d2 = 2*(6.0+log((double)size)/dlog2);
	if(d2 > 53.0){
	    cout <<"erfc: error; contact author"<<endl;
	    exit(0);
	}

	//Allocation of needed memory
	ePow = new double*[size];
	if(ePow){
	    for(int i=0; i<size; i++){
		ePow[i] = new double[mp::mp5];
		if(ePow[i] == NULL){
		    cout << "QuadErf::erfc: Not enough Memory for e^() table" << endl;
		    exit(-1);
		}else ePow[i][0] = mp::mp5;
	    }
	}else
	{
	    cout << "QuadErf::erfc: Not enough Memory for e^() table" << endl;
	    exit(-1);
	}

	//t1 = - alpha ** 2
	alpha2 = alpha*alpha;
	c_mpdmc(-alpha2, sk0);

	//t2 = exp (t1)
	c_mpexp(sk0,sk1);

	// t3 = t2 ** 2
	c_mpmul(sk1,sk1,sk2);

	//t4 = 1.d0
	c_mpeq(one,sk3);

	for(int i=0; i< size; i++){
	    //t4 = t2 * t4
	    c_mpmul(sk1,sk3,sk4);
	    c_mpeq(sk4,sk3);
	    
	    //ePow[i] = t4
	    c_mpeq(sk3, ePow[i]);

	    //t2 = t2 * t3  (oddSequence)
	    c_mpmul(sk1,sk2,sk4);
	    c_mpeq(sk4,sk1);

	}
    }

    //t1 = 0.d0
    c_mpdmc(0.0,sk0);

    //t2 = t ** 2
    c_mpmul(t,t, sk1);

    //t3 = exp (-t2)
    c_mpneg_q(sk1,sk4);
    c_mpexp(sk4,sk2);
    
    int ic;
    for(int k=0; k<size; k++){
	//t5 = ePow[k] / (k2 ** 2 * alpha ** 2 + t2)
	int k2 = k+1;
	c_mpdmc((k2*k2)*alpha2, sk4);
	c_mpadd(sk4,sk1,sk5);
	c_mpdiv(ePow[k], sk5, sk4);

	//t1 = t1 + t5
	c_mpadd(sk0,sk4,sk5);
	c_mpeq(sk5,sk0);

	//if (abs (t5) < t4) break;
	c_mpabs(sk4,sk5);
	c_mpltt(sk5,eps,&ic);
	if(ic >0) break;
    }
    

    //erfc = t3 * alpha * t / PI * (1.d0 / t2 + 2.d0 * t1) &
    //       + 2.d0 / (1.d0 - exp (2.d0 * PI * t / alpha))
    c_mpmul_qd(sk2,alpha,sk3);
    c_mpmul(sk3,t,sk4);
    c_mpdiv(sk4,mp_real::_pi.mpr, sk3);
    c_mpdiv(one,sk1,sk4);
    c_mpmul_qd(sk0,2.0,sk5);
    c_mpadd(sk4,sk5,sk6);
    c_mpmul(sk3,sk6,sk4);
    c_mpmul_qd(mp_real::_pi.mpr,2.0,sk3);
    c_mpmul(sk3,t,sk5);
    c_mpdiv_qd(sk5,alpha,sk6);
    c_mpexp(sk6,sk3);
    c_mpsub(one,sk3,sk5);
    c_mpdmc(2.0,sk6);
    c_mpdiv(sk6,sk5,sk3);

    double* ans = new double[mp::mp5];
    if(ans){
	ans[0] = mp::mp5;
	c_mpadd(sk4,sk3,ans);
    }else{
	cout << "QuadErf::erfc: Not enough Memory create answer" << endl;
	exit(-1);
    }
    erfMPIDErrorHandle(arg);
    delete sk0;
    delete sk1;
    delete sk2;
    delete sk3;
    delete sk4;
    delete sk5;
    return mp_real(ans);
}
*/