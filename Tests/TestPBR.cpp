#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std ;

// Sources
#include "Point/Point.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

// gMLQMC
#include <gMLQMC/gMLQMC.hpp>
#include <gMLQMC/samplers/uniform.hpp>
#include <gMLQMC/samplers/IPL.hpp>
#include <gMLQMC/samplers/halton.hpp>
#include <gMLQMC/strategy/full.hpp>

// Fourier Representation
#include "BoundaryRepresentations/Fourier.hpp"

int main(int argc, char* argv[]) {
	int NumberofParameters = 10 ;
	ParametricBoundaryRepresentation PBR(NumberofParameters,R_Fourier,R_der_Fourier) ;
	vector <double> y ;
	y.assign(NumberofParameters,0.0) ; 
	PBR.SetParameters(y) ;
	// Quadrature
	int NumberofPoints = 24 ;
	int NumberofCycles =  1 ; 
	PBR.SetQuadrature(NumberofPoints,NumberofCycles) ;
	Point P ;
	PBR.GetBoundaryRepresentation(P,0.0) ;
	cout << P[0] << " " << P[1] << "\n" ;
	cout << PBR.ComputeBoundaryMeasure() << "\n" ;
	typedef ParametricBoundaryRepresentation PPP;
	return 0 ;
}

