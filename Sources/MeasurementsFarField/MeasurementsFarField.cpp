#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>
#include <map> 

// Sources
#include "QoIBoundary.hpp"
#include "Point/Point.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

using namespace std ;

MeasurementsFarField::MeasurementsFarField() {
	
}

MeasurementsFarField::MeasurementsFarField(BoundaryPotentials & BP_, int NumberOfPoints_) {
	BP = &BP_ ;
	NumberOfPoints = NumberOfPoints_ ;
	Xcoord.assign(NumberOfPoints,0.0) ;
	Ycoord.assign(NumberOfPoints,0.0) ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		ArcParam[i] = (double)(i)/(double)(NumberOfControlPoints) ;
	}
}

void QoIBoundary::ComputePositions(ParametricBoundaryRepresentation & PBR) {
	Point P ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		PBR.GetBoundaryRepresentation(P,ArcParam[i]) ;
		PointsBoundary[2*i  ] = P[0] ;
		PointsBoundary[2*i+1] = P[1] ;
	}
}

QoIBoundary & QoIBoundary::operator=(const QoIBoundary & QoI) {
	NumberOfControlPoints = QoI.NumberOfControlPoints ;
	ArcParam = QoI.ArcParam ;  
	PointsBoundary = QoI.PointsBoundary ;
	return *this;	
}

QoIBoundary QoIBoundary::operator+(const QoIBoundary & QoI) const {
	QoIBoundary NewQoI(QoI.NumberOfControlPoints) ;
	NewQoI.ArcParam = QoI.ArcParam ;
	assert(NumberOfControlPoints == QoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		NewQoI.PointsBoundary[i] = PointsBoundary.at(i) + (QoI.PointsBoundary).at(i) ;	
	}
	return NewQoI;
}	

QoIBoundary QoIBoundary::operator-(const QoIBoundary & QoI) const {
	QoIBoundary NewQoI(QoI.NumberOfControlPoints) ;
	NewQoI.ArcParam = QoI.ArcParam ;
	assert(NumberOfControlPoints == QoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		(NewQoI.PointsBoundary).at(i) = PointsBoundary.at(i) - (QoI.PointsBoundary).at(i) ;	
	}
	return NewQoI;
}

QoIBoundary QoIBoundary::operator*(double a) const {
	QoIBoundary NewQoI(NumberOfControlPoints) ; 
	assert(NumberOfControlPoints == QoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		(NewQoI.PointsBoundary).at(i) = ((NewQoI.PointsBoundary).at(i)) * a ;	
	}
	return NewQoI ;
}

void QoIBoundary::PrintData() {
	for(int i = 0; i<NumberOfControlPoints; i++) {
		cout << PointsBoundary[2*i  ] << "    " << PointsBoundary[2*i+1] << "\n" ;
	}
	cout << "\n" ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		cout << ArcParam[i] << "\n" ;
	}
	
}

QoIBoundary operator*(double a, QoIBoundary const & QoI) {
	return QoI*a ;
}
	