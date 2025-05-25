//Implenetation of functions in util.h
#include <cmath>
#include <iostream>
#include <arprec/c_mp.h>  // should clean this up
#include "util.h"

using std::abs;
using std::cout;
using std::endl;

/*
 *  This function returns true when there is an error_no error
 *  If the errorStr value is set it will be printed in stdout
 *     Value of error will also be printed out along with this message
 *  error value will also set in error
 */
bool isThereErrorMPIER(int & error, const char* errorStr)
{
    error = mp::error_no;
    //c_mpgetpar(6/*error_no*/, &error,0);
    if(error > 0)
    {
	if (errorStr)
	    cout << errorStr << "\nError Value = " << error <<endl;
	return true;
    }
    return false;
}


/*
 *  This function returns true when there is an error_no error
 *  If the errorStr value is set it will be printed in stdout
 */
bool isThereErrorMPIER(const char* errorStr)
{
    int error;
    return isThereErrorMPIER(error, errorStr);
}


/*
 * returns the DP approximation to log10 (a)
 *
 */
double dplog10q(mp_real a)
{
    static const double log10_2  = log10(2.0);

    int ia;
    double da;
    c_mpmdc(a.mpr,&da,&ia);
    if(da == 0.0)
	return -9999.0;
    
    return log10(da < 0 ? -da : da) + ia*log10_2;
}


/*
 * For input MP value a, this routine returns DP b and integer ib such that 
 * a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.
 */
void decmdq(const mp_real& a, double& b, int& ib)
{
    static const double xlt = 0.3010299956639812;
    double da;
    int ia;
    c_mpmdc(a.mpr,&da, &ia);

    if(da != 0.0){
	double t1 = xlt*ia + log10(abs(da));
	ib = (int)t1;
	if(t1 < 0.0) ib -= 1;
	
	//b = sign(pow(10.0,t1-ib), da);
	b = pow(10.0, t1-ib);
	if(da < 0)
	    b = -b;

	return;
    }

    b = 0.0;
    ib = 0;
}

