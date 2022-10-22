#ifndef SPECTRALQUANTITY_HEADERDEF 
#define SPECTRALQUANTITY_HEADERDEF

#include <cmath> 
#include <complex>
#include <map>
#include <vector>

using namespace std ;

class SpectralQuantity {

private:
	int nModes ;
	int Samples ;
	map < int , complex <double> > SQ ;

public:
	SpectralQuantity() ;
	SpectralQuantity(int nModes_, int Samples_ = 0) ;
	SpectralQuantity(const SpectralQuantity & SQ_) ;
	void ConvertToSpectral(complex < double > (*f)(double t, void * Input), void * Input_) ;
	complex <double> & operator[](int i) ;
	complex <double> operator()(double t) ; 
	complex <double> GetScaled(double t, complex <double> ScalingFactor) ;
	SpectralQuantity & operator=(const SpectralQuantity & S) ;
	SpectralQuantity operator+(const SpectralQuantity & S) const ;
	SpectralQuantity operator-(const SpectralQuantity & S) const ;
	SpectralQuantity operator*(double a) const ;
	int GetNumberOfModes() ;
	friend SpectralQuantity operator*(double a, SpectralQuantity const & S) ;
} ;	
	SpectralQuantity operator*(double a, SpectralQuantity const & S) ;
#endif