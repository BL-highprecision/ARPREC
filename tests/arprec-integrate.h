/*
 * ArprecIntegrate.h
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
#ifndef _LBL_G0V_ARPREC_INTEGRATE_H_
#define _LBL_G0V_ARPREC_INTEGRATE_H_

#include "integrate.h"

template <typename Space>
class ArprecIntegrate: public Integrate<Space>
{
 public:
    /*
     * reInitializes the instance, if there is a need to do so
     * including the child class
     */
    virtual void reInit(int phases, long int AbWtSize, int neps1,
			int neps2, int debug);
    
    //gets the phases
    inline int getPhases(){return _phases;}

    //gets the precWord1
    inline int getPrecWd1(){return _precWd1;}

    //gets the precWord2
    inline int getPrecWd2(){return _precWd2;}


    /*
     *  Following are virtual functions intented to 
     *  get values from child classes
     */
    
    //Gets the Max Phases Child class is configured with
    virtual int getMaxPhases()=0;
    
    //Gets the Max Precision in words for primary precision
    //Child class is configured with
    virtual int getMaxPrecWd1()=0;

    //Gets the Max Precision in words for secondary precision
    //Child class is configured with
    virtual int getMaxPrecWd2()=0;
    

    /*
     * sets the Phases to the passed value
     * if false is returned reinit must be called
     * to set the phases to the desired value
     */
    bool setPhases(int phases);

    /*
     * sets the precision in words for primary
     * Returns false if the value is bigger than
     * Max prec
     */
    bool setPrecWd1(int prec);

    /*
     * sets the precision in words for secondary
     * Returns false if the value is bigger than
     * Max prec
     */
    bool setPrecWd2(int prec);

 protected:
    ArprecIntegrate(int phases, int neps1, int neps2,
		    int debug, int nerror=0, int ierror=0);

 
    //Resizes the Abscissas and weight arrays
    virtual void reSizeAbWt(long int AbWtSize) = 0;
    
    int _phases;  //Keeps track of the number of phases
    int _neps1;   //Log10 of the primary tolerance
    int _precWd1; //precesion in words, corresponding to primary tolerence
    int _neps2;   //Log10 of the secondary tolerance
    int _precWd2; //precesion in words, corresponding to secondary tolerence

    int _ndebug;  //used for printing debug info
    //Error tracking variables
    int _nerror;
    int _ierror;
};

#include "arprec-integrate.cpp"

#endif
