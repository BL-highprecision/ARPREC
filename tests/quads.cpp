/*
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 7, 2004
 * Modified:      July 2, 2004
 * Version:       1.0
 */
#include <cstdio>
#include <ctime>
#include <arprec/mp_real.h>
#include "quad-erf.h"
#include "quad-gs.h"
#include "quad-ts.h"
#include "util.h"

using std::cout;
using std::endl;

int nwords1=1, nwords2=1;

/*
 * List of function used for testing the quadrature
 *
 */

mp_real fun01(const mp_real &t)
{
    return t*log(1.0+t);
}

mp_real fun02(const mp_real &t)
{
    return sqr(t)*atan(t);
}

mp_real fun03(const mp_real &t)
{
    return exp(t)*cos(t);
}

mp_real fun04(const mp_real &t)
{
    mp_real t2 = sqr(t);
    mp_real temp = sqrt(2.0 + t2);
    return atan(temp)/((1+t2)*temp);
}

mp_real fun05(const mp_real &t)
{
    return sqrt(t)*log(t);
}

mp_real fun06(const mp_real &t)
{
    return sqrt(1.0-sqr(t));
}

mp_real fun07(const mp_real &t)
{
    mp::mpsetprecwords(nwords2);
    mp_real t1(1.0-t);
    mp::mpsetprecwords(nwords1);
    return t/sqrt(t1*(1.0+t));
}

mp_real fun08(const mp_real &t)
{
    return sqr(log(t));
}

mp_real fun09(const mp_real &t)
{
    mp_real t1 = cos(t);
    if (t1 > 0.0)
	return log(t1);
    return mp_real(0.0);
}

mp_real fun10(const mp_real &t)
{
//    cout << "fun10 arg =\n" << t;
    if(t > .25*mp_real::_pi){
	mp::mpsetprecwords(nwords2);
	mp_real t1= .5*mp_real::_pi -t;
	mp::mpsetprecwords(nwords1);
//	cout << "fun10 transArg =\n" << t1;
	return 1.0/sqrt(tan(t1));
    }
    return sqrt(tan(t));
}

mp_real fun11(const mp_real &t)
{
    return 1.0/(1.0 -2.0*t + 2.0*sqr(t));
}

mp_real fun12(const mp_real &t)
{
    mp::mpsetprecwords(nwords2);
    mp_real t1(1.0-t);
    mp::mpsetprecwords(nwords1);
    return exp(1.0-1.0/t)/sqrt(pow(t,3)*t1);
}

mp_real fun13(const mp_real &t)
{
    return exp(-.5*sqr(1.0/t-1.0))/sqr(t);
}

mp_real fun14(const mp_real &t)
{
    mp_real t1(1.0/t-1.0);
    return exp(-1*t1)*cos(t1)/sqr(t);
}

mp_real fun15(const mp_real &t)
{
    return sin(t)/t;
}

mp_real fun16(const mp_real &t)
{
    return pow(t,7)*sin(1/t);
}


typedef mp_real (* testFunc)(const mp_real &ty) ;

testFunc funArray[]={fun01,fun02,fun03,fun04,fun05,fun06,fun07,fun08,
		     fun09,fun10,fun11,fun12,fun13,fun14,fun15,fun16};

enum QUADS {QUAD_GS=1, QUAD_ERFC=2, QUAD_TS=3};


/*
 *  USAGE OF MAIN
 *  all parameter to main will be nubmers
 *  First number will tell which Quadrature routine to use,
 *     so must be between 1-3
 *          1 => QuadGS
 *          2 => QuadErfc
 *          3 => QuadTS
 *  Second number will set the Quadrature levels(phases)
 *  Third Will tell the size of the Table
 *  Fourth will be the number of digits in primary precision
 *  Fifth will be the number of digits in secondary precision
 *  Sixth will be the debug level setting
 */
int main(int argc, char **argv){
#define NUM_PARAMETRS 6
    int argt[NUM_PARAMETRS] ={QUAD_TS, 10, 10, 400, 800, 2};
    argc--;
    if(argc > 0){
	if(argc > NUM_PARAMETRS)
	    argc = NUM_PARAMETRS;
	for(int i=0; i< argc; i++)
	    argt[i] = atoi(argv[i+1]);
    }

    QUADS quad = (QUADS)argt[0];
    int phases=argt[1];
    long int AbWtSize=argt[2];
    int neps1=-argt[3];
    int neps2=-argt[4];
    int debug=argt[5];

    mp::mp_init(800);
    mp::mpsetoutputprec(400);
    time_t t1 = time(NULL),t2;

    ArprecIntegrate<mp_real> *integrator;

    switch(quad)
    {
	case QUAD_GS:
	    if(argc<=1)
		phases--;
	    if(argc<=4)
		neps2 = neps1;
	    integrator = new QuadGS(phases, AbWtSize, neps1, neps2, debug);
	    break;

	case QUAD_ERFC:
	    integrator = new QuadErf(phases, AbWtSize, neps1, neps2, debug);
	    break;

	case QUAD_TS: default:
	    integrator = new QuadTS(phases, AbWtSize, neps1, neps2, debug);
    }
    
    t2 = time(NULL);

    cout <<"Quadrature initialization completed: cpu time ="<< t2-t1<<endl;
    cout << "Continuous functions on finite itervals:" << endl;
    
    nwords1 = integrator->getPrecWd1();
    nwords2 = integrator->getPrecWd2();

    mp::mpsetprecwords(nwords2);
    mp_real x1(0.0), x2(1.0);
    mp::mpsetprecwords(nwords1);
    mp_real val(0.0),a_val(0.0);
    double d1; int n1;


    cout << "Problem 1: Int_0^1 t*log(1+t) dt = 1/4" <<endl;
    t1=time(NULL);
    val = integrator->integrate(funArray[0],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = .25;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 2: Int_0^1 t^2*arctan(t) dt = (pi - 2 + 2*log(2))/12" <<endl;
    t1=time(NULL);
    val = integrator->integrate(funArray[1],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = (mp_real::_pi - 2.0 + 2*mp_real::_log2)/12.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 3: Int_0^(pi/2) e^t*cos(t) dt = 1/2*(e^(pi/2) - 1)" <<endl;
    x1 = 0.0;
    x2 = .5*mp_real::_pi;
    t1=time(NULL);
    val = integrator->integrate(funArray[2],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = .5*(exp(.5*mp_real::_pi) - 1.0);
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 4: Int_0^1 arctan(sqrt(2+t^2))/((1+t^2)sqrt(2+t^2)) dt = 5*Pi^2/96" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[3],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = 5.0*mp_real::_pi*mp_real::_pi/96.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Continuous functions on finite itervals, but non-diff at an endpoint"
	 << endl;

    cout << "Problem 5: Int_0^1 sqrt(t)*log(t) dt = -4/9" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[4],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = mp_real(-4.0)/9.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;

    cout << "Problem 6: Int_0^1 sqrt(1-t^2) dt = pi/4" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[5],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = mp_real::_pi/4.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Functions on finite intervals with integrable singularity at an endpoint.\n"
	 << "Problem 7: Int_0^1 t/sqrt(1-t^2) dt = 1" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[6],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = 1.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 8: Int_0^1 log(t)^2 dt = 2" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[7],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = 2.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 9: Int_0^(pi/2) log(cos(t)) dt = -pi*log(2)/2" <<endl;
    x1 = 0.0;
    x2 = .5*mp_real::_pi;
    t1=time(NULL);
    val = integrator->integrate(funArray[8],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = -.5*mp_real::_pi*log(mp_real(2.0));
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 10: Int_0^(pi/2) sqrt(tan(t)) dt = pi*sqrt(2)/2" <<endl;
    mp::mpsetprecwords(nwords2);
    x1 = 0.0;
    x2 = .5*mp_real::_pi;
    mp::mpsetprecwords(nwords1);
    t1=time(NULL);
    val = integrator->integrate(funArray[9],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = mp_real::_pi*sqrt(mp_real(2.0))/2.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Functions on an infinite interval (requiring a two-step solution\n"
	 << "Problem 11: Int_0^inf 1/(1+t^2) dt = pi/2"<<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[10],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = mp_real::_pi/2.0;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 12: Int_0^inf e^(-t)/sqrt(t) dt = sqrt(pi)" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[11],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = sqrt(mp_real::_pi);
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 13: Int_0^inf e^(-t^2/2) dt = sqrt(pi/2)" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[12],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = sqrt(mp_real::_pi/2.0);
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "\n\nOscillatory functions on an infinite interval.\n"
	 <<  "Problem 14: Int_0^inf e^(-t)*cos(t) dt = 1/2" <<endl;
    x1 = 0.0;
    x2 = 1.0;
    t1=time(NULL);
    val = integrator->integrate(funArray[13],x1,x2);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = .5;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;


    cout << "Problem 15: Int_0^inf sin(t)/t = pi/2" <<endl;
    x1 = 0.0;
    x2 = mp_real::_pi;
    t1=time(NULL);
    val = integrator->integrate(funArray[14],x1,x2);
    x2 = 1.0/mp_real::_pi;
    val += 40320.0*integrator->integrate(funArray[15],x1,x2)
	-1.0/mp_real::_pi + 2.0/pow(mp_real::_pi,3)
	-24.0/pow(mp_real::_pi,5) + 720.0/pow(mp_real::_pi,7);
    t2=time(NULL);
    cout << "Quadrature Completed: CPU time =" << t2-t1 
	 << "\nResult=\n" << val;
    a_val = .5*mp_real::_pi;
    decmdq(a_val - val, d1, n1);
    cout << "Actual error = " << d1 << "x10^"<<n1<<"\n"<<endl;
    cout <<"Prob 15 error may be 40,000 X higher than estimated error"<<endl;

#undef NUM_PARAMETRS
    delete integrator;
    mp::mp_finalize();
    return 0;
}

