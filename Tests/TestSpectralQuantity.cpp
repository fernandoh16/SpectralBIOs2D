#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std ;

// Sources
#include "Point/Point.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

//
#include "BoundaryRepresentations/Circle.hpp"

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		// cout << P[0] << " " << P[1] << "\n" ;
		return -(PBR->NormOfTangentVector(t)) * E * exp(complex<double>(0.0,WaveNumber * (Direction & P))) ;
	}
} ;

complex<double> GetPlaneWave(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->Wave(t) ;
}

int main(int argc, char* argv[]) {
	int NumberOfModes = 10 ;
	int nSamples = 200 ;
	double WaveNumber = 1.0 ;
	double Magnitude = 1.0 ;
	double DirectionX = 1.0 ;
	double DirectionY = 0.0 ;
	int NumberOfPoints = 24 ;
	int NumberOfCycles = 1 ;
	int NumberOfParameters = 100 ;
	// === Parametric Boundary Representation ===
	ParametricBoundaryRepresentation PBR(NumberOfParameters,R_Fourier,R_der_Fourier) ;
	vector <double> y ;
	y.assign(NumberOfParameters,1.0) ; 
	PBR.SetParameters(y) ;
	// Quadrature 
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	// === RHS ===
	SpectralQuantity RHS(NumberOfModes,nSamples) ;
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
	int NumberOfControlPoints = 50 ;
	double EvalPoint = 0.0 ;
	cout << GetPlaneWave(EvalPoint,&PW) << "\n" ;
	/*
	for(int i = -NumberOfModes; i<=NumberOfModes; i++) {
		cout << RHS[i].real() << "  " ;
		cout << RHS[i].imag() << "\n" ;
	}
	*/
	for(int i = 0; i<NumberOfControlPoints; i++) {
		EvalPoint = (double)(i) / (double)(NumberOfControlPoints) ;
		cout << GetPlaneWave(EvalPoint,&PW) << "  " << RHS.GetScaled(EvalPoint,1.0) << "\n" ;
	}
	return 0 ;
}

