#ifndef QoIBOUNDARY_HPP
#define QoIBOUNDARY_HPP

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

// Sources
#include "Point/Point.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

using namespace std ;

class QoIBoundary {
private:
	int      NumberOfControlPoints ;
	vector <double>       ArcParam ;
	vector <double> PointsBoundary ;
	
public:
	QoIBoundary() ;
	QoIBoundary(ParametricBoundaryRepresentation & PBR_, int NumberOfControlPoints_) ;
	void ComputePositions(ParametricBoundaryRepresentation & PBR_, Point & P) ;
	/*
	QoIBoundary & operator=(const QoIBoundary & QoI) ;
	QoIBoundary operator+(const QoIBoundary & QoI) const ;
	QoIBoundary operator-(const QoIBoundary & QoI) const ;
	QoIBoundary operator*(double a) const ;
	void PrintData() ;
	friend QoIBoundary operator*(double a, QoIBoundary const & QoI) ;
	*/
} ;
	// QoIBoundary operator*(double a, QoIBoundary const & QoI) ;
#endif

				