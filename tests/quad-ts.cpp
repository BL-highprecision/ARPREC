/*
 * quad-ts.cpp
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

//Implementation of QuadTS
#include <cmath>
#include <iostream>
#include "quad-ts.h"
#include "util.h"

using std::cout;
using std::endl;

//Static Variable initialization
int      QuadTS::_phasesMax    = 0;
long int QuadTS::_AbWtSize     = 0;
long int QuadTS::_AbWtMax      = 0;
long int QuadTS::_quadTSCount  = 0;
mp_real* QuadTS::_abscissas    = NULL;
mp_real* QuadTS::_weight       = NULL;
int      QuadTS::_precWord1    = 0;
int      QuadTS::_precWord2    = 0;
int      QuadTS::_sneps1       = 0;
int      QuadTS::_sneps2       = 0;


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
void QuadTS::reSizeAbWt(long int AbWtSize)
{
    if(_phasesMax < _phases)
	_phasesMax = _phases;

    if(AbWtSize < 12*2)
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
	    cout << "QuadTS::reSizeAbWt: Not enough memory" <<endl;
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
	cout <<"QuadTS::reSizeAbWt: Tanh-sinh quadrature initialization"
	     <<endl;

    /*
     *  Some inits before actual computation of tables begin
     *
     */
    mp::mpsetprecwords(_precWord2);
    mp_real eps2 = pow(mp_real(10.0), _sneps2);
    mp::mpsetprecwords(_precWord1);
    mp_real p2(mp_real::_pi/2),
	t1,t2,t3,t4, u1,u2;
    double h = std::ldexp(1.0, -_phasesMax);
#define IPRINT 1000

    for(int k=_AbWtSize; k < _AbWtMax; k++)
    {
	//   xk(k) = 1 - tanh (u1) = 1 /(e^u1 * cosh (u1))
	//   wk(k) = u2 / cosh (u1)^2
	//   where u1 = pi/2 * cosh (t1), u2 = pi/2 * sinh (t1)
	t1 = k*h;
	t2 = exp(t1);
	u1 = .5*p2*(t2 + 1.0/t2);
	u2 = .5*p2*(t2 - 1.0/t2);
	t3 = exp(u2);
	t4 = .5*(t3 + 1.0/t3);
	_abscissas[k] = 1.0/(t3*t4);
	_weight[k]    = u1/(t4*t4);

	if(_ndebug >= 2 && k%IPRINT == 0)
	{
	    cout<<"\t\t"<<k<<"\t"<<_AbWtMax<<endl;
	    //cout<<"abscissas[i]=\n"<<_abscissas[k]
	    //<<"weight[i]=\n"<<_weight[k];
	}

	 if(isThereErrorMPIER(_ierror))
	 {
	     _nerror = _ierror + 100;
	     cout <<"reSizeAbWt: Error in quadratore initialization; code"
		  << _nerror <<endl;
	     return;
	 }
	 if(_weight[k] < eps2)
	 {
	     _AbWtSize = k;
	     if(_ndebug >= 2)
		 cout <<"reSizeAbWt: Tabale space used = " 
		      << _AbWtSize << endl;
	     return;
	 }
    }
    cout<< "reSizeAbWt: Table space parameter is too small; value ="
	<< _AbWtMax <<endl;
    _AbWtSize = _AbWtMax;
}

mp_real QuadTS::integrate(mp_real func(const mp_real &x),
			  const mp_real &x1, const mp_real &x2)
{
    static const int izx = 5;
    int outprec = mp::mpgetoutputprec();
    mp::mpsetoutputprec(56);

    mp::mpsetprecwords(_precWd2);
    mp_real ax = .5*(x2-x1);
    mp_real bx = .5*(x2+x1);

    mp::mpsetprecwords(_precWd1);
    mp_real sum(0.0),s1(0.0),s2(0.0),s3(0.0),
	fm1(0.0),fm2(0.0), err(0.0), h(1.0),
	xx1, xx2, eps1, eps2,
	t1, t2, tw1, tw2;
    const mp_real c10(10.0);
    const mp_real ceps1(pow(c10, _neps1));
    const mp_real ceps2(pow(c10, _neps2));
    
    int phases = _phases+1, k1, k2, iz1,iz2;
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
	k1  = ip[_phases-k];
	k2  = ip[_phases-k+1];

	iz1 = iz2 = 0;

	//Evaluate function at level k in x, avoiding unnecessary computation.
	for(long int i = 0; i < _AbWtSize; i += k1)
	{
	    if(i%k2 != 0 || k == 1)
	    {
		mp::mpsetprecwords(_precWord2);
		mp_real xt1(1.0 - _abscissas[i]);
		mp_real xx1(bx - ax*xt1);
		mp_real xx2(bx + ax*xt1);
		bool log1 = xx1 > x1 && iz1 < izx;
		bool log2 = xx2 < x2 && iz2 < izx;
		mp::mpsetprecwords(_precWord1);

		if(!log1 && !log2)
		    break;

		if(log1)
		{
		    t1 = func(xx1);
		    if(isThereErrorMPIER(_ierror) || _nerror > 0)
		    {
			if(_ierror > 0) _nerror = 100 + _ierror;
			cout << "QuadTS::integrate: Error in quadrature"
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
			cout << "QuadTS::integrate: Error in quadrature"
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
	    cout << "QuadTS::integrate: Iteration=" << k << " of " << _phases
		 << "; est error = 10^" << nint(dble(dplog10q(abs(err))))
		 << "; approx value = " << s1 << endl;

	if(k > 3 && err < eps1){
	    mp::mpsetoutputprec(outprec);
	    delete ip;
	    return s1;
	}
	
	if(k >= 3 && err < eps2)
	{
	    cout << "QuadTS::integrate: Estimated error = 10^" 
		 << nint(dble(dplog10q(abs(err)))) << "\n"
		 << "Increase secondary prec (Ndigits2) for greater accuracy. "
		 << "Current value =" << -_neps2 << endl;
	    mp::mpsetoutputprec(outprec);
	    delete ip;
	    return s1;
	}
    }
    
    cout << "QuadTS::integrate: Estimated error = 10^"
	 << nint(dble(dplog10q(abs(err)))) << "\n"
	 << "Increae QuadLevel for greater accuracy. "
	 << "Current value =" << _phases << endl;

    mp::mpsetoutputprec(outprec);
    delete ip;
    return s1;
}