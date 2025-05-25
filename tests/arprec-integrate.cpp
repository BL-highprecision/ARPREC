/*
 * ArprecIntegrate.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 25, 2004
 * Modified:      June 25, 2004
 * Version:       1.0
 */
//Implementation of ArprecIntegrate


#ifndef _LBL_G0V_ARPREC_INTEGRATE_CC_
#define _LBL_G0V_ARPREC_INTEGRATE_CC_

#include "arprec-integrate.h"
#include <arprec/mp.h>

/*
 *  Constructor to do the initialization
 *
 */
template <typename Space>
ArprecIntegrate<Space>::ArprecIntegrate(int phases, int neps1, int neps2,
					int debug, int nerror, int ierror):
    _phases(phases), _neps1(neps1), _neps2(neps2),
    _ndebug(debug), _nerror(nerror), _ierror(ierror)
{
    mp::mpsetprec(-_neps2);
    _precWd2  = mp::mpgetprecwords();
    mp::mpsetprec(-_neps1);
    _precWd1 = mp::mpgetprecwords();
}

/*
 * reInitializes the instance, if there is a need to do so
 * including the child class
 */
template <typename Space>
void ArprecIntegrate<Space>::reInit(int phases, long int AbWtSize, int neps1,
				    int neps2, int debug)
{
    _phases = phases;

    if(_neps2 != neps2)
    {
	_neps2 = neps2;
	mp::mpsetprec(-_neps2);
	_precWd2  = mp::mpgetprecwords();
    }
    if(_neps1 != neps1)
    {
	_neps1 = neps1;
	mp::mpsetprec(-_neps1);
	_precWd1  = mp::mpgetprecwords();
    }
    _ndebug  = debug;
    _nerror = 0;
    _ierror = 0;

    reSizeAbWt(AbWtSize);
}


/*
 * sets the Phases to the passed value
 * if false is returned reinit must be called
 * to set the phases to the desired value
 */
template <typename Space>
inline bool ArprecIntegrate<Space>::setPhases(int phases)
{
    if(phases <= getMaxPhases() && phases > 0)
    {
	_phases = phases;
	return true;
    }
    return false;
}


/*
 * sets the precision in words for primary
 * Returns false if the value is bigger than
 * Max prec
 */
template <typename Space>
bool ArprecIntegrate<Space>::setPrecWd1(int prec)
{
    if(prec <= getMaxPrecWd1() && prec >0)
    {
	_precWd1 = prec;
	mp::mpsetprecwords(prec);
	_neps1   = -mp::mpgetprec();
	return true;
    }
    return false;
}


/*
 * sets the precision in words for secondary
 * Returns false if the value is bigger than
 * Max prec
 */
template <typename Space>
bool ArprecIntegrate<Space>::setPrecWd2(int prec)
{
    if(prec <= getMaxPrecWd2() && prec >0)
    {
	_precWd2 = prec;
	mp::mpsetprecwords(prec);
	_neps2   = -mp::mpgetprec();
	return true;
    }
    return false;
}

#endif
