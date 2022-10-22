#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std ;

#include "Point/Point.hpp"
#include "Quadrature/Quadrature.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"

void R(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] = cos(2.0 * M_PI * (double)(Index) * t) ;
		// P[1] = sin(2.0 * M_PI * (double)(Index) * t) ;
		P[0] = 10.0 * cos(2.0 * M_PI * t) ;
		P[1] = sin(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}

void R_der(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] =-sin(2.0 * M_PI * (double)(Index) * t) ;
		// P[1] = cos(2.0 * M_PI * (double)(Index) * t) ;
		P[0] =-2.0 * 10.0 * M_PI * sin(2.0 * M_PI * t) ;
		P[1] = 2.0 * M_PI * cos(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}

int main(int argc, char* argv[]) {
	int NumberofParameters = 10 ;
	ParametricBoundaryRepresentation PBR(NumberofParameters,R,R_der) ;
	vector <double> y ;
	y.assign(NumberofParameters,0.0) ; 
	PBR.SetParameters(y) ;
	/*
	// Quadrature
	int NumberofPoints = 24 ;
	int NumberofCycles =  1 ; 
	PBR.SetQuadrature(NumberofPoints,NumberofCycles) ;
	Point P ;
	PBR.GetBoundaryRepresentation(P,0.0) ;
	// cout << P[0] << " " << P[1] << "\n" ;
	// cout << PBR.ComputeBoundaryMeasure() << "\n" ;
	// Spectral BIOs
	SpectralBIOs S(2,"Helmholtz",1.0,&PBR) ;
	//S.BuildSpectralBIOsOld() ; 
	S.BuildSpectralBIOs(1500) ; 
	//S.ComputeDifference() ;
	*/
	return 0 ;
}

