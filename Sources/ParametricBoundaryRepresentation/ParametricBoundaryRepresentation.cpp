#include <cmath> 
#include <complex>
#include <iostream>
#include <cassert>
#include <map> 

using namespace std;

#include "ParametricBoundaryRepresentation.hpp"
#include "Point/Point.hpp"
#include "Quadrature/Quadrature.hpp"

struct BoundaryMeasureStruct {
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	double BoundaryTransformationNorm(double t) {
		(*PBR).GetBoundaryRepresentationDerivative(P,t) ;
		return length(P) ;
	}
} ;

double BoundaryMeasure(double t, void * Input) {
	BoundaryMeasureStruct * BMS = (BoundaryMeasureStruct *) Input ;
	return (BMS->BoundaryTransformationNorm(t)) ;
}

ParametricBoundaryRepresentation::ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_)) {
	NumberofParameters = NumberofParameters_ ;
	R = R_ ;
	yParameters.assign(NumberofParameters,0.0) ; 
}

ParametricBoundaryRepresentation::ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_), void (*R_der_)(Point & P_, double t_, int Index_)) {
	NumberofParameters = NumberofParameters_ ;
	R = R_ ;
	R_der = R_der_ ;
	yParameters.assign(NumberofParameters,0.0) ; 
}

ParametricBoundaryRepresentation::ParametricBoundaryRepresentation(int NumberofParameters_, void (*R_)(Point & P_, double t_, int Index_), void (*R_der_)(Point & P_, double t_, int Index_), void (*R_sec_der_)(Point & P_, double t_, int Index_)) {
	NumberofParameters = NumberofParameters_ ;
	R = R_ ;
	R_der = R_der_ ;
	R_sec_der = R_sec_der_ ;
	yParameters.assign(NumberofParameters,0.0) ; 
}

void ParametricBoundaryRepresentation::SetQuadrature(int NumberofPoints_, int NumberofCycles_) {
	Q.SetQuadrature("1D",NumberofPoints_,NumberofCycles_) ;
}

int ParametricBoundaryRepresentation::GetNumberofParameters() {
	return NumberofParameters ;
}

void ParametricBoundaryRepresentation::SetParameters(vector <double> & yParameters_) {
	assert(yParameters.size() == yParameters.size()) ;
	yParameters = yParameters_ ;
}

void ParametricBoundaryRepresentation::GetBoundaryRepresentation(Point & P_, double t) {
	Point P ;
	R(P,t,0) ;
	P_ = P ;
	for (int j = 0; j<NumberofParameters; j++) {
		R(P,t,j+1) ;
		P_ = P_ + (P * (2.0 * yParameters[j] - 1.0)) ;
	}
}

void ParametricBoundaryRepresentation::GetBoundaryRepresentationDerivative(Point & P_, double t) {
	Point P ;
	R_der(P,t,0) ;
	P_ = P ;
	for (int j = 0; j<NumberofParameters; j++) {
		R_der(P,t,j+1) ;
		P_ = P_ + (P * (2.0 * yParameters[j] - 1.0)) ;
	}
}

void ParametricBoundaryRepresentation::GetBoundaryRepresentationSecondDerivative(Point & P_, double t) {
	Point P ;
	R_sec_der(P,t,0) ;
	P_ = P ;
	for (int j = 0; j<NumberofParameters; j++) {
		R_sec_der(P,t,j+1) ;
		P_ = P_ + (P * (2.0 * yParameters[j] - 1.0)) ;
	}
}

double ParametricBoundaryRepresentation::ComputeBoundaryMeasure() {
	BoundaryMeasureStruct BM ;
	BM.PBR  = this ;
	return Q.ComputeQuadrature(BoundaryMeasure,(void*)&BM,0.0,1.0) ;
}

double ParametricBoundaryRepresentation::NormOfTangentVector(double t) {
	Point P ;
	GetBoundaryRepresentationDerivative(P,t) ;
	return length(P) ;
}

void ParametricBoundaryRepresentation::GetNormalVector(Point & P_, double t) {
	Point P(0.0,0.0) ;
	GetBoundaryRepresentationDerivative(P,t) ;
	P_[0] = +P[1]/length(P) ;
	P_[1] = -P[0]/length(P) ;
}