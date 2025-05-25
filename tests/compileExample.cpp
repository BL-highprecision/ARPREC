

// simple exaple of ARPREC usage to illustrate linking process
// Alex Kaiser, LBNL, 6/3/2010



#define numDigitsPrecision1 500


#include <iostream>
#include <arprec/mp_real.h>
#include <arprec/mp_int.h>

using namespace std; 

int main() {
	
	// initialization should be set to desired precision plus two
	mp::mp_init(numDigitsPrecision1 + 2); 
	mp::mpsetprec(numDigitsPrecision1 ); 
	mp::mpsetoutputprec(numDigitsPrecision1 ); 
	cout.precision(numDigitsPrecision1 ) ; 

	// simple read example
	/*
	mp_real readTest ; 
	cin >> readTest ; 	
	cout << "readTest = " << readTest << endl ;
	 */ 
	
	// simple demo
	mp_real x = "1.0" ;  
	x /= 3.0 ;
	mp_real y ; 
	y = pow( mp_real(2.0) , 3 ) ; 
	cout << "y = " << y << endl;
	cout << "x = " << x << endl; 
	
	
	mp_real a ; 
	mp_real b = mp_real("0.1");
	
	a = sqrt(b);
	cout << " sqrt(0.1) = " << a << endl;
	cout << " sqrt(0.1) * sqrt(0.1) = " << a * a << endl; 
	
	
	mp::mp_finalize();
	return 0;
}


