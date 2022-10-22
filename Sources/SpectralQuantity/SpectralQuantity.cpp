#include <cmath> 
#include <complex>
#include <iostream>
#include <map> 
#include <cassert> 

// FFT
#include <fftw3.h>

// Sources
#include "SpectralQuantity.hpp"

using namespace std;

SpectralQuantity::SpectralQuantity() {
	nModes  = 0  ;
	Samples = 0 ;
}

SpectralQuantity::SpectralQuantity(int nModes_, int Samples_) {
	nModes  = nModes_  ;
	Samples = Samples_ ;
	for (int i = -nModes; i<=nModes; i++) {
		SQ[i] = complex<double>(0.0,0.0) ;
	}
	if(Samples_ == 0) {
		Samples_ = nModes ;
	}
}

SpectralQuantity::SpectralQuantity(const SpectralQuantity & SQ_) {
	nModes  = SQ_.nModes ;
	Samples = SQ_.Samples ;
	SQ = SQ_.SQ ;
}

void SpectralQuantity::ConvertToSpectral(complex < double > (*f)(double t, void * Input), void * Input_) {
	// assert((Samples>=nModes)) ;
	// assert(Samples % 2 = 1) ;
	int NewSamples = 2 * Samples + 1 ;
	// FFT Objects
  	fftw_complex * Input ;
  	fftw_complex * Output ;
  	//
  	Input  = new fftw_complex[NewSamples] ;
  	Output = new fftw_complex[NewSamples] ;
  	fftw_plan P ;
  	//
  	double S ;
  	complex<double> V ;
  	for (int i = 0; i<NewSamples; i++) {
  		S = (double)(i)/(double)(NewSamples) ;
  		V = (*f)(S,Input_) ;
    	Input[i][0] = V.real() ;
    	Input[i][1] = V.imag() ;
  	}
	// 
  	P = fftw_plan_dft_1d(NewSamples,Input,Output,FFTW_FORWARD,FFTW_ESTIMATE) ;
  	fftw_execute(P) ;
  	fftw_destroy_plan(P) ;
  	/*
    for(int i = 0; i<NewSamples; i++) {
  		cout << i<< "   "<<Output[i][0] / NewSamples << "   " << Output[i][1] / NewSamples << "\n" ;
  	}
  	cout << "\n\n" ;
  	*/
  	// Store data in SQ
  	int Index ;
  	for (int i = 0; i<=Samples; i++) {
  		if(i == 0) {
  			SQ[0] = 1.0/(double)(NewSamples) * complex <double>(Output[0][0],Output[0][1]) ;
  		}
  		else {
			Index  = 2 * Samples - i + 1;
			SQ[-i] = 1.0/(double)(NewSamples) * complex <double>(Output[Index][0],Output[Index][1]) ;
  			Index  = i ;
			SQ[+i] = 1.0/(double)(NewSamples) * complex <double>(Output[Index][0],Output[Index][1]) ;
		}
  	}
  	//
  	/*
  	for(int i = -Samples; i<=Samples; i++) {
  		cout <<i<< "   " << SQ[i] << "\n" ;
  	}
  	*/
    fftw_free(Input) ;
  	fftw_free(Output) ;
}

complex<double> & SpectralQuantity::operator[](int i) {
	assert((i<=nModes) && (i>=-nModes)) ;
	return SQ[i] ;
}

complex<double> SpectralQuantity::operator()(double t) {
	complex<double> Value  = complex<double>(0.0) ;
	for(int i = -nModes; i<=nModes; i++) {
		Value = Value + SQ[i] * exp(complex<double>(0.0,2.0 * M_PI * (double)(i) * t)) ;
	}
	return Value ;
}

complex<double> SpectralQuantity::GetScaled(double t, complex <double> ScalingFactor) {
	complex<double> Value  = complex<double>(0.0) ;
	for(int i = -nModes; i<=nModes; i++) {
		Value = Value + SQ[i] * exp(complex<double>(0.0,2.0 * M_PI * (double)(i) * t)) ;
	}
	return (Value / ScalingFactor) ;
}

SpectralQuantity & SpectralQuantity::operator=(const SpectralQuantity & S) {
	nModes  =  S.nModes ;
	Samples = S.Samples ;  
	SQ = S.SQ ;
	return *this;	
} 

SpectralQuantity SpectralQuantity::operator+(const SpectralQuantity & S) const {
	SpectralQuantity NewSQ(S.nModes,S.Samples) ; 
	for(int i = -nModes; i<=nModes; i++) {
		NewSQ[i] = SQ.at(i) + (S.SQ).at(i) ;	
	}
	return NewSQ;
}

SpectralQuantity SpectralQuantity::operator-(const SpectralQuantity & S) const {
	SpectralQuantity NewSQ(S.nModes,S.Samples) ; 
	for(int i = -nModes; i<=nModes; i++) {
		NewSQ[i] = SQ.at(i) - (S.SQ).at(i) ;	
	}
	return NewSQ ;
}

SpectralQuantity SpectralQuantity::operator*(double a) const {
	SpectralQuantity NewSQ(nModes,Samples) ; 
	for(int i = -nModes; i<=nModes; i++) {
		NewSQ[i] = SQ.at(i) * a ;	
	}
	return NewSQ;
}

int SpectralQuantity::GetNumberOfModes() {
	return nModes ;
}

SpectralQuantity operator*(double a, SpectralQuantity const & S) {
	return S*a ;
}


