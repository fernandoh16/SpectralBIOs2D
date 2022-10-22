#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>
#include <map> 

// Sources
#include "QoISolution.hpp"
#include "Point/Point.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"

using namespace std ;

QoISolution::QoISolution() {
	
}

QoISolution::QoISolution(int NumberOfControlPoints_) {
	NumberOfControlPoints = NumberOfControlPoints_ ;
	ArcParam.assign(NumberOfControlPoints,0.0) ;
	Solution.assign(NumberOfControlPoints,complex<double>(0.0,0.0)) ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		ArcParam[i] = (double)(i)/(double)(NumberOfControlPoints) ;
	}
}

void QoISolution::ComputeSolution(SpectralQuantity & SQ) {
	Point P ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		Solution[i] = SQ(ArcParam[i]) ;
	}
}

QoISolution & QoISolution::operator=(const QoISolution & QoI) {
	NumberOfControlPoints = QoI.NumberOfControlPoints ;
	ArcParam = QoI.ArcParam ;  
	Solution = QoI.Solution ;
	return *this;	
}

QoISolution QoISolution::operator+(const QoISolution & QoI) const {
	QoISolution NewQoI(QoI.NumberOfControlPoints) ;
	NewQoI.ArcParam = QoI.ArcParam ;
	assert(NumberOfControlPoints == QoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		NewQoI.Solution[i] = Solution.at(i) + (QoI.Solution).at(i) ;	
	}
	return NewQoI;
}	

QoISolution QoISolution::operator-(const QoISolution & QoI) const {
	QoISolution NewQoI(QoI.NumberOfControlPoints) ;
	NewQoI.ArcParam = QoI.ArcParam ;
	assert(NumberOfControlPoints == QoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		(NewQoI.Solution).at(i) = Solution.at(i) - (QoI.Solution).at(i) ;	
	}
	return NewQoI;
}

QoISolution QoISolution::operator*(double a) const {
	QoISolution NewQoI(NumberOfControlPoints) ; 
	assert(NumberOfControlPoints == NewQoI.NumberOfControlPoints) ; 
	for(int i = 0; i<2*NumberOfControlPoints; i++) {
		(NewQoI.Solution).at(i) = ((NewQoI.Solution).at(i)) * a ;	
	}
	return NewQoI ;
}

void QoISolution::PrintData() {
	for(int i = 0; i<NumberOfControlPoints; i++) {
		cout << Solution[i] << "\n" ;
	}
	cout << "\n" ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		cout << ArcParam[i] << "\n" ;
	}
	
}

QoISolution operator*(double a, QoISolution const & QoI) {
	return QoI*a ;
}
	