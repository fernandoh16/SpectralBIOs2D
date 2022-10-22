#ifndef PARAMETRICREPRESENTATION_HEADERDEF 
#define PARAMETRICREPRESENTATION_HEADERDEF

#include <cmath> 
#include <complex>
#include <map>
#include <vector>

using namespace std ;

#include "Point/Point.hpp"
#include "Quadrature/Quadrature.hpp"

class ParametricBoundaryRepresentation {

private:
	int NumberofParameters ;
	vector <double> yParameters ;
	void (*R)(Point & P_, double t, int Index) ;
	void (*R_der)(Point & P, double t, int Index) ;
	void (*R_sec_der)(Point & P, double t, int Index) ;
	Quadrature Q ; 

public:
	ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_)) ;
	ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_), void (*R_der_)(Point & P_, double t_, int Index_)) ;
	ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_), void (*R_der_)(Point & P_, double t_, int Index_), void (*R_sec_der_)(Point & P_, double t_, int Index_)) ;
	void SetQuadrature(int  NumberofPoints_, int NumberofCycles_) ;
	int GetNumberofParameters() ;
	void SetParameters(vector <double> & yParameters_) ;
	void GetBoundaryRepresentation(Point & P_, double t) ;
	void GetBoundaryRepresentationDerivative(Point & P_, double t) ;
	void GetBoundaryRepresentationSecondDerivative(Point & P_, double t) ;
	double ComputeBoundaryMeasure() ;
	double NormOfTangentVector(double t) ;
	void GetNormalVector(Point & P_, double t) ;
};
	
#endif