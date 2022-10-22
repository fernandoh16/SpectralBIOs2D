#define _USE_MATH_DEFINES
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
		// P[0] = 10.0 * cos(2.0 * M_PI * t) ;
		// P[1] = sin(2.0 * M_PI * t) ;
		P[0] = cos(2.0 * M_PI * t) + 0.65 * cos(4.0 * M_PI * t) - 0.65 ;
		P[1] = 1.5 * sin(2.0 * M_PI * t) ;
	}
	else {
		// P[0] = 0.0 ;
		// P[1] = 0.0 ;
	}
}

void R_der(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] =-10.0 * 2.0 * M_PI * sin(2.0 * M_PI * t) ;
		// P[1] = 2.0 * M_PI * cos(2.0 * M_PI * t) ;
		P[0] = -2.0 * M_PI * sin(2.0 *M_PI*t) - 4.0 * M_PI * 0.65 * sin(2.0 * M_PI * t) ;
		P[1] = 1.5 * 2.0 * M_PI * cos(2.0 * M_PI * t) ;
	}
	else {
		// P[0] = 0.0 ;
		// P[1] = 0.0 ;
	}
}

void R_sec_der(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] =-10.0 * pow(2.0 * M_PI,2.0) * cos(2.0 * M_PI * t) ;
		// P[1] =-pow(2.0 * M_PI,2.0) * sin(2.0 * M_PI * t) ;
		P[0] = -pow(2.0 * M_PI,2.0) * cos(2.0 *M_PI*t) - pow(4.0 * M_PI,2.0) * 0.65 * cos(2.0 * M_PI * t) ;
		P[1] = -1.5 * pow(2.0 * M_PI,2.0) * sin(2.0 * M_PI * t) ;
	}
	else {
		// P[0] = 0.0 ;
		// P[1] = 0.0 ;
	}
}

int main(int argc, char* argv[]) {
	int NumberofParameters = 10 ;
	int NumberOfModes = 10;
	int nSamples = 500;
	string Problem = "Helmholtz";
	double WaveNumber = 1.0;
	ParametricBoundaryRepresentation PBR(NumberofParameters,R,R_der,R_sec_der) ;
	vector <double> y ;
	y.assign(NumberofParameters,0.0) ; 
	PBR.SetParameters(y) ;
	//
	SpectralBIOs S(NumberOfModes,&PBR,nSamples,Problem,WaveNumber) ; 
	S.BuildSpectralV() ;
	S.BuildSpectralK() ;
	S.BuildMassMatrix() ;
	return 0 ;
}

