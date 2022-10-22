#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>
#include <map> 

// Sources
#include "BoundaryPotentials.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "Quadrature/Quadrature.hpp"

// Boost
#include <boost/math/special_functions/hankel.hpp> 

using namespace std ;

struct BoundaryPotentialsStruct {
	Point P ;
	Point P_Boundary ;
	Point Nu ;
	ParametricBoundaryRepresentation * PBR ;
	double WaveNumber ;
	complex<double> SingleLayer(double t) {	
		PBR->GetBoundaryRepresentation(P_Boundary, t) ;
		return complex<double>(0.0,1.0)/4.0 * boost::math::cyl_hankel_1(0, WaveNumber * length(P,P_Boundary)) ;
	}
	complex<double> DoubleLayer(double t) {	
		PBR->GetBoundaryRepresentation(P_Boundary, t) ;
		PBR->GetBoundaryRepresentationDerivative(Nu,t) ;
		return -complex<double>(0.0,WaveNumber * ((P-P_Boundary) & Nu)/(length(P,P_Boundary) * length(Nu))/4.0) * boost::math::cyl_hankel_1(1, WaveNumber * length(P,P_Boundary)) ;
	}
	complex<double> FarFieldSingleLayer(double t) {	
		PBR->GetBoundaryRepresentation(P_Boundary, t) ;
		// Check !!
		return exp(-complex<double>(0.0,M_PI/4.0)) / sqrt(8.0 * M_PI * WaveNumber) * exp(complex<double>(0.0,-WaveNumber * (P & P_Boundary)/length(P))) ;
	}
	complex<double> FarFieldDoubleLayer(double t) {	
		// === TODO ===
		return complex<double>(0.0,0.0) ;
	}
} ;

complex<double> SL(double t, void * Input) {
	BoundaryPotentialsStruct * BP = (BoundaryPotentialsStruct *) Input ;
	return BP->SingleLayer(t) ;
}

complex<double> DL(double t, void * Input) {
	BoundaryPotentialsStruct * BP = (BoundaryPotentialsStruct *) Input ;
	return BP->DoubleLayer(t) ;
}

complex<double> FFSL(double t, void * Input) {
	BoundaryPotentialsStruct * BP = (BoundaryPotentialsStruct *) Input ;
	return BP->FarFieldSingleLayer(t) ;
}

complex<double> FFDL(double t, void * Input) {
	BoundaryPotentialsStruct * BP = (BoundaryPotentialsStruct *) Input ;
	return BP->FarFieldDoubleLayer(t) ;
}

BoundaryPotentials::BoundaryPotentials(SpectralQuantity & SQ_,  ParametricBoundaryRepresentation & PBR_, int nSamples_ , double WaveNumber_) {
	SQ = &SQ_ ;
	PBR = &PBR_ ;
	nSamples = nSamples_ ; 
	WaveNumber = WaveNumber_ ;
}

complex <double> BoundaryPotentials::SingleLayerPotential(Point & P) {
	BoundaryPotentialsStruct BPS ;
	BPS.P = P ;
	BPS.PBR = PBR ;
	BPS.WaveNumber = WaveNumber ;
	//
	SpectralQuantity SQ_SLP(SQ->GetNumberOfModes(),nSamples) ;
	SQ_SLP.ConvertToSpectral(SL,&BPS) ;
	int Index = min(SQ->GetNumberOfModes(),(nSamples-1)/2) ;
	//
	complex<double> Value = complex<double>(0.0,0.0) ;
	for(int i = -Index; i<=Index; i++) {
		Value = Value + SQ_SLP[i] * (*SQ)[i] ;
	}
	return Value ;
}

complex <double> BoundaryPotentials::DoubleLayerPotential(Point & P) {
	BoundaryPotentialsStruct BPS ;
	BPS.P = P ;
	BPS.PBR = PBR ;
	BPS.WaveNumber = WaveNumber ;
	//
	SpectralQuantity SQ_DLP(SQ->GetNumberOfModes(),nSamples) ;
	SQ_DLP.ConvertToSpectral(DL,&BPS) ;
	int Index = min(SQ->GetNumberOfModes(),(nSamples-1)/2) ;
	//
	complex<double> Value = complex<double>(0.0,0.0) ;
	for(int i = -Index; i<=Index; i++) {
		Value = Value + SQ_DLP[i] * (*SQ)[i] ;
	}
	return Value ;
}

complex <double> BoundaryPotentials::FarFieldSingleLayerPotential(Point & P) {
	BoundaryPotentialsStruct BPS ;
	BPS.P = P ;
	BPS.PBR = PBR ;
	BPS.WaveNumber = WaveNumber ;
	//
	SpectralQuantity SQ_FFSLP(SQ->GetNumberOfModes(),nSamples) ;
	SQ_FFSLP.ConvertToSpectral(FFSL,&BPS) ;
	//
	complex<double> Value = complex<double>(0.0,0.0) ;
	int Index = SQ->GetNumberOfModes() ;
	for(int i = -Index; i<=Index; i++) {
		Value = Value + SQ_FFSLP[i] * (*SQ)[-i] ;
	}
	return Value ;
}

complex <double> BoundaryPotentials::FarFieldDoubleLayerPotential(Point & P) {
	BoundaryPotentialsStruct BPS ;
	BPS.P = P ;
	BPS.PBR = PBR ;
	BPS.WaveNumber = WaveNumber ;
	//
	SpectralQuantity SQ_FFDLP(SQ->GetNumberOfModes(),nSamples) ;
	SQ_FFDLP.ConvertToSpectral(FFDL,&BPS) ;
	int Index = min(SQ->GetNumberOfModes(),(nSamples-1)/2) ;
	//
	complex<double> Value = complex<double>(0.0,0.0) ;
	for(int i = -Index; i<=Index; i++) {
		Value = Value + SQ_FFDLP[i] * (*SQ)[i] ;
	}
	return Value ;
}

				