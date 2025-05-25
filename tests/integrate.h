/*
 * Integrate.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2004
 *
 *
 * Creation Date: June 16, 2004
 * Modified:      June 24, 2004
 * Version:       1.0
 */

#ifndef _LBL_GOV_INTEGRATE_H_
#define _LBL_GOV_INTEGRATE_H_

/*
 *  Parent class for integration qudratures
 *
 */
template <typename Space> 
class Integrate
{
 public:
    static const double pi;// = 3.141592653589793238; //pi in double precision

    virtual Space integrate(Space func(const Space &x), const Space &x1,
			    const Space &x2)=0;
    virtual ~Integrate() { }
};

/*
 * Static const init as it cannot be handled in the class def by
 * some of the c++ compilers
 */

//pi in double precision
template<typename Space>
const double Integrate<Space>::pi = 3.141592653589793238; 

#endif
