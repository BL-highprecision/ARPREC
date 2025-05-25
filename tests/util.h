/*
 * Util.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 10, 2004
 * Modified:      June 24, 2004
 * Version:       1.0
 */

/*
 * This contains some common functions used by the Quadrature Classes
 */

#ifndef _LBL_GOV_UTIL_H_
#define _LBL_GOV_UTIL_H_

#ifndef NULL
#define NULL 0
#endif

#include <arprec/mp_real.h>

/*
 *  This function returns true when there is an error_no error
 *  If the errorStr value is set it will be printed in stdout
 *  error value will also set in error
 */
bool isThereErrorMPIER(int & error, const char* errorStr = NULL);

/*
 *  This function returns true when there is an error_no error
 *  If the errorStr value is set it will be printed in stdout
 */
bool isThereErrorMPIER( const char* errorStr = NULL);


/*
 *  calculates the DP approximation to Log10(a)
 *
 */
double dplog10q(mp_real a);

/*
 * For input MP value a, this routine returns DP b and integer ib such that 
 * a = b * 10^ib, with 1 <= abs (b) < 10 for nonzero a.
 */
void decmdq(const mp_real& a, double& b, int& ib);

#endif
