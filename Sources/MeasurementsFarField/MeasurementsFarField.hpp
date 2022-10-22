#ifndef MEASUREMENTSFARFIELD_HPP
#define MEASUREMENTSFARFIELD_HPP

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

// Sources
#include "Point/Point.hpp"
#include "BoundaryPotentials/BoundaryPotentials.hpp"

using namespace std ;

class MeasurementsFarField {
private:
	BoundaryPotentials * BP ;
	int NumberOfPoints ;
	vector < double > Xcoord ;
	vector < double > YCoord ;
	vector < complex < double > > M ;
	
public:
	MeasurementsFarField() ;
	MeasurementsFarField(BoundaryPotentials & BP, NumberOfPoints) ;
	MeasurementsFarField & operator=(const QoIBoundary & QoI) ;
	MeasurementsFarField operator+(const QoIBoundary & QoI) const ;
	MeasurementsFarField operator-(const QoIBoundary & QoI) const ;
	MeasurementsFarField operator*(double a) const ;
	// void PrintData() ;
	// void ReadData() ;
	friend MeasurementsFarField operator*(double a, MeasurementsFarField const & QoI) ;
} ;
	MeasurementsFarField operator*(double a, MeasurementsFarField const & QoI) ;
#endif

				