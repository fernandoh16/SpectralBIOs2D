#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

// Sources
#include "Point/Point.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "QoIBoundary/QoIBoundary.hpp"

// Fourier Representation
#include "BoundaryRepresentations/Fourier.hpp"

using namespace std ;

int main(int argc, char* argv[]) {
	// === PBR ===
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
	//cout << P[0] << " " << P[1] << "\n" ;
	// cout << PBR.ComputeBoundaryMeasure() << "\n" ;
	// === QoIBoundary ===
	int NumberOfControlPoints = 10 ;
	QoIBoundary QoI_1(NumberOfControlPoints) ;
	QoIBoundary QoI_2 ;
	QoIBoundary QoI_3 ;
	//
	QoI_1.ComputePositions(PBR) ;
	QoI_1.PrintData() ;
	QoI_2 = QoI_1 ;
	QoI_2.PrintData() ;
	QoI_3 = QoI_1 - QoI_2 ; 
	QoI_3.PrintData() ;
	return 0 ;
}

