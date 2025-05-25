/*
 * quad-gs.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 7, 2004
 * Modified:      June 26, 2004
 * Version:       1.0
 */

//Implementation of QuadGS
#include <cmath>
#include <iostream>
#include "quad-gs.h"
#include "util.h"

using std::cout;
using std::endl;

//Static Variable initialization
int      QuadGS::_phasesMax    = 0;
int      QuadGS::_AbWtSize     = 0;
int      QuadGS::_AbWtMax      = 0;
int      QuadGS::_quadGSCount  = 0;
mp_real* QuadGS::_abscissas    = NULL;
mp_real* QuadGS::_weight       = NULL;
int      QuadGS::_precWord1    = 0;
int      QuadGS::_precWord2    = 0;
int      QuadGS::_sneps1       = 0;
int      QuadGS::_sneps2       = 0;


/*
 *  ABSTRACT FUNCTION
 *  Resizes the Abscissas and weight arrays
 *
 *  This subroutine initializes the quadrature arays xk and wk for Gaussian
 *  quadrature.  It employs a Newton iteration scheme with a dynamic precision
 *  level.  The argument nq2, which is the space allocated for wk and xk in
 *  the calling program, should be at least 8 * 2^nq1 + 100, although a higher
 *  value may be required, depending on precision level.  Monitor the space
 *  figure given in the message below during initialization to be certain.
 *  2004-06-25
 */
void QuadGS::reSizeAbWt(long int AbWtSize)
{
    //Variable used in minimizing the calculation
    //if not all the things has to be recomputed
    static int phasesTrav = 1, n=3*2*2, nby2=3*2, j=1;

    if(_phasesMax < _phases)
	_phasesMax = _phases;

    if(AbWtSize < 12*4 + 100)
    {
	AbWtSize = 12;
	//AbWtSize = 8;
	for(int i =1; i<_phasesMax; i++)
	    AbWtSize *= 2;
	AbWtSize += 100;
    }
    //Need to recompute everything as precesion has been increased
    if(_precWord1 < _precWd1){
	_AbWtSize = _AbWtMax=  0;
	phasesTrav = 1;
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
	    cout << "QuadGS::reSizeAbWt: Not enough memory" <<endl;
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
	cout <<"QuadGS::reSizeAbWt: Gaussian quadrature initialization"
	     <<endl;

    /*
     *  Some inits before actual computation of tables begin
     *
     */
    int ik0 = 100, i;
    mp_real eps(std::ldexp(1.0, -96));
    mp_real r, t1, t2, t3, t4, t5;

    if(phasesTrav == 1){
	_weight[0] = _abscissas[0] = 0.0;
	_weight[1] = _phasesMax;
	_abscissas[1] = ik0;
	i = ik0;

	for(int j=2; j <= ik0; j++)
	    _abscissas[j] = _weight[j] = 0.0;
    }else i = std::max(ik0,_AbWtSize-1);

    //phasesTrav has been set to its value already
    for(; phasesTrav <= _phasesMax; )
    {
	if(_ndebug >= 2)
	    cout<<"\t\t"<<phasesTrav<<"\t"<<i<<"\t"<<_AbWtMax<<endl;

	for(;j<=nby2; j++)
	{
 	    if(++i >= _AbWtMax){
		cout<< "reSizeAbWt: Table space parameter is too small; "
		    << "value =" << _AbWtMax <<endl;
		_AbWtSize = _AbWtMax;
		return;
	    }

	    // Compute a double precision estimate of the root
	    int nwp = 3;
	    mp::mpsetprecwords(nwp);
	    r = cos( (pi * (j-.25) ) / (n+.5) );

	    //Compute the j-th root of the n-degree 
	    //Legendre polynomial using Newton's iteration.
	    for(;;)//Infinite Loop
	    {
		t1 = 1.0;
		t2 = 0.0;
		for(int j1 = 1; j1 <= n; j1++)
		{
		    t3 = t2;
		    t2 = t1;
		    t1 = ((2*j1-1)*r*t2-(j1-1)*t3)/j1;
		}
		t4 = n*(r*t1 - t2)/(r*r - 1.0);
		t5 = r;
		r -= t1/t4;

		if(nwp==3)
		{
		    mp_real temp(r - t5);
		    if(abs(temp) <= eps)
		    {
			nwp *= 2;
			nwp--;
			nwp = std::min(nwp,_precWord1);
			mp::mpsetprecwords(nwp);
		    }
		}else if(nwp < _precWord1)
		{
		    nwp *= 2;
		    nwp--;
		    nwp = std::min(nwp,_precWord1);
		    mp::mpsetprecwords(nwp);
		}else break;
	    }
	    
	    _abscissas[i] = r;
	    t4 = n*(r*t1-t2) /( r*r - 1.0);
	    _weight[i] = 2.0/((1.0-r*r)*t4*t4);
	    if(isThereErrorMPIER(_ierror))
	    {
		_nerror = _ierror + 100;
		cout <<"reSizeAbWt: Error in quadratore initialization; code="
		     << _nerror <<endl;
		nby2 = 3*2;
		n    = 3*2*2;
		j    = 1;
		phasesTrav = 1;
		return;
	    }
	}
	nby2 = n;
	n *= 2;
	j = 1;
	_abscissas[++phasesTrav] = i;
    }
    nby2 = 3*2;
    n    = 3*2*2;
    _AbWtSize = i+1;
    cout <<"reSizeAbWt: Table space used = " << _AbWtSize << endl;

}

mp_real QuadGS::integrate(mp_real func(const mp_real &x),
			  const mp_real &x1, const mp_real &x2)
{
    int outprec = mp::mpgetoutputprec();
    mp::mpsetoutputprec(56);

    mp::mpsetprecwords(_precWd1);

    mp_real ax = .5*(x2-x1);
    mp_real bx = .5*(x2+x1);
    mp_real sum(0.0),s1(0.0),s2(0.0),s3(0.0),
	fm1(0.0),fm2(0.0), err(0.0),
	xx1, xx2, eps1, eps2, eps1x,
	t1, t2, t3;
    const mp_real c10(10.0), zero(0.0);
    const mp_real ceps1(pow(c10, _neps1));
    const mp_real ceps2(pow(c10, _neps2));
    
    int n = 3*2,i, nerr; 

    for(int k=1; k<=_phases; k++)
    {
	int nby2 = n;
	n *= 2;
	s3 = s2;
	s2 = s1;

	fm1 = 0.0;
	fm2 = 0.0;
	sum = 0.0;
	i = (int)dble(_abscissas[k]);

	for(int j=1; j <= nby2; j++)
	{
	    i++;
	    xx1 = bx -ax*_abscissas[i];
	    xx2 = bx +ax*_abscissas[i];
	    
	    if(xx1 > x1)
	    {
		t1 = func(xx1);
		if(isThereErrorMPIER(_ierror) || _nerror > 0)
		{
		    if(_ierror > 0)_nerror = 100 + _ierror;
		    cout << "QuadGS::integrate: Error in quadrature"
			"calculation; code =" << _nerror << endl;
		    mp::mpsetoutputprec(outprec);
		    return mp_real(0.0);
		}
	    }else t1 = 0.0;

	    if(xx2 < x2 && j+k > 2)
	    {
		t2 = func(xx2);
		if(isThereErrorMPIER(_ierror) || _nerror > 0)
		{
		    if(_ierror > 0) _nerror = 100 + _ierror;
		    cout << "QuadGS::integrate: Error in quadrature"
			 << "calculation; code =" << _nerror << endl;
		    mp::mpsetoutputprec(outprec);
		    return mp_real(0.0);
		}
	    }else t2 = 0.0;
	    
	    t3 =_weight[i] * (t1 + t2);
	    sum += t3;

	    t3 = abs(t3);
	    t2 = abs(t2);
	    t1 = abs(t1);
	    fm1 = std::max(fm1, t3);
	    fm2 = std::max(fm2, t2);
	    fm2 = std::max(fm2, t1);
	}

	//Intergral = s1 AND Error estimation calculation follows
	s1 = ax*sum;
	eps1  = fm1*ceps1;
	eps1x = fm1* ceps1; 
	eps2  = fm2*ceps2;
	double d1 = dplog10q(abs(s1 - s2));
	double d2 = dplog10q(abs(s1 - s3));
	double d3 = dplog10q(eps1) - 1;
	double d4 = dplog10q(eps2) - 1;


	if(k<=2)
	{
	    nerr = 0;
	    err = 1.0;
	}else if(d1 == -9999.0)
	{
	    nerr = -9999;
	    err  = 0.0;
	}else
	{
	    double val1 = d1*d1/d2;
	    double val2 = 2*d1;
	    double max = std::max(val1,val2);
	    max  = std::max(max,d3);
	    max  = std::max(max,d4);
	    nerr = (int)nint(std::min(0.0,max));
	    err  = pow(c10, nerr);
	}

	if(_ndebug >= 2)
	    cout << "QuadGS::integrate: Iteration=" << k << " of " << _phases
		 << "; est error = 10^" << nint(dble(dplog10q(abs(err))))
		 << "; approx value = " << s1 << endl;

	if(k > 3 && err < eps1){
	    mp::mpsetoutputprec(outprec);
	    return s1;
	}

	if(k >= 3 && err < eps2)
	{
	    cout << "QuadGS::integrate: Estimated error = 10^" 
		 << nint(dble(dplog10q(abs(err)))) << "\n"
		 << "Increase secondary prec (Ndigits2) for greater accuracy. "
		 << "Current value =" << -_neps2 << endl;
	    mp::mpsetoutputprec(outprec);
	    return s1;
	}

    }

    cout << "QuadGS::integrate: Estimated error = 10^"
	 << nint(dble(dplog10q(abs(err)))) << "\n"
	 << "Increae QuadLevel for greater accuracy. "
	 << "Current value =" << _phases << endl;

    if( err > 1.0e-20)
	cout <<"quadgs: Poor results may be due to singularities at "
	     <<"endpoints.\nIf so, try the erf or tanh-sinh "
	     <<"quadrature routines (Quadtype = 2 or 3)."
	     <<endl;

    mp::mpsetoutputprec(outprec);
    return s1;
}