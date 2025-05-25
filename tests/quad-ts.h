/*
 * QuadTS.h
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
class QuadTS;

#ifndef _LBL_GOV_QUADTS_H_
#define _LBL_GOV_QuADTS_H_

#include "arprec-integrate.h"
#include <arprec/mp_real.h>

/*
 *  This is the implementation of the quadrature routine 'quaderf', which employs
 *  the error function.  The function quaderf is suitable to integrate
 *  a function that is continuous, infinitely differentiable and integrable on a
 *  finite open interval.  It can also be used for certain integrals on
 *  infinite intervals, by making a suitable change of variable -- see below.
 *  While this routine is not quite as efficient as quadgs for functions that 
 *  are regular on a closed interval, it can be used for functions with an
 *  integrable singularity at one or both of the endpoints.
 *
 *  The function(s) to be integrated is(are) passed in as a parameter(function 
 *  pointers) along with the domain of intergraion
 *
 *  Note that an integral of a function on an infinite interval can be
 *  converted to an integral on a finite interval by means of a suitable
 *  change of variable.  Example (here the notation "inf" means infinity):
 *
 *  Int_0^inf f(t) dt  =  Int_0^1 f(t) dt + Int_1^inf f(t) dt
 *                     =  Int_0^1 f(t) dt + Int_0^1 f(1/t)/t^2 dt
 */
class QuadTS: public ArprecIntegrate<mp_real>
{
 public:
    /*
     *  neps* are the values that control the precision of the answer
     *  neps1  Log10 of the primary tolerance.
     *  neps2  Log10 of the secondary tolerance.
     *
     *  phases is the max number of phases in the quadrature routine.
     *  adding one  increases (possibly doubles) the number of accurate 
     *  digits in the result, but also roughly doubles the run time.
     *  phases must be atlest 2, if phases < 2 it will be set to 2.
     *  Size of the Abscissas and weight array is set to the default
     *  value of 12*2^phases
     *
     *  AbWtSize is the size of the Abscissus and weight array
     */
    QuadTS(int phases=2, long int AbWtSize=48, int neps1=-400,
	   int neps2=-800, int debug=2): 
	ArprecIntegrate<mp_real>(phases,neps1,neps2,debug)
    {
	_quadTSCount++;  //Used to decide whether or not to delete the arrays
	reSizeAbWt(AbWtSize);
    }

    /*
     *  Destructor for the Class
     *  it is used mainly for deleting the abscissas and weight
     *  array. It will delete if there were no quadTS class
     *  There will be no QuadTS class if quadTSCount = 0
     */
    virtual ~QuadTS()
    {
	_quadTSCount--;
	//Deleting the static arrays if this is the last instance 
	//of QuadTS to be deleted
	if(_quadTSCount == 0)
	{
	    if (_abscissas) delete [] _abscissas;
	    if(_weight) delete [] _weight;
	    _abscissas = _weight = NULL;
	    _AbWtSize = _AbWtMax = _precWord1 = _precWord2 = 0;
	}
    }

    
    /*
     * Abstract Functions
     *
     */

    /*
     * Returns the value of integrating func in interval
     * [x1, x2]
     */
    mp_real integrate(mp_real func(const mp_real &x), 
		      const mp_real &x1, const mp_real &x2);

    //gets the max phases
    inline int getMaxPhases(){return _phasesMax;}

    //gets the Actual precWord1
    //Precision of the Abscissas and weight
    inline int getMaxPrecWd1(){return _precWord1;}

    //gets the Actual precWord2
    //Precision of the Abscissas and weight
    inline int getMaxPrecWd2(){return _precWord2;}

    //gets the Abscissus array
    inline static const mp_real* getAbscissus(){return _abscissas;}

    //gets the weight array
    inline static const mp_real* getWeight(){return _weight;}

    //gets the size of the Abscissus and weight array
    inline static long int getAbWtSize(){return _AbWtSize;}

 protected:
    //ABSTRACT FUNCTION
    //Resizes the Abscissas and weight arrays
    void reSizeAbWt(long int AbWtSize);

 private:

    //Static variables as they wil be same for all QuadErf classes
    //These variables will also be used to determine when content
    //of the static arrays will be changed
    static int _phasesMax;         //The number of phases
    static long int _AbWtSize;     //Number of values in the arrays
    static long int _AbWtMax;      //size of the Abscissas and weight Array
    static mp_real* _abscissas;    //Pointer to Abscissas array
    static mp_real* _weight;       //Pointer to Weight array

    static int _precWord1;         //The word precesion with which the table
    static int _precWord2;         //was computed

    static int _sneps1;            //Smallest Log 10 primary and secondary
    static int _sneps2;            //tolerence QuadErf initiated with
    
    static long int _quadTSCount; //Used to decide when to delete the arrays

};

#endif
