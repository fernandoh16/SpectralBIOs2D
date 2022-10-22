#ifndef QOI_HPP
#define QOI_HPP
#pragma once

using namespace std ;

#include "Point/Point.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

struct QuantityOfInterest {
	ParametricBoundaryRepresentation & PBR ;
	int      NumberOfControlPoints ;
	vector <double>       ArcParam ;
	QuantityOfInterest(ParametricBoundaryRepresentation & PBR_, int NumberOfControlPoints_): 
					PBR(PBR_) {
		NumberOfControlPoints = NumberOfControlPoints_ ;
		ArcParam.assign(NumberOfControlPoints,0.0) ;
		for(int i = 0; i<NumberOfControlPoints; i++) {
			ArcParam[i] = (double)(i)/(double)(NumberOfControlPoints) ;
		}
	}
	vector< double > operator()(vector < double > & y) {
		PBR.SetParameters(y) ;
		Point P ;
		vector <double> PointsBoundary ;
		PointsBoundary.assign(2 * NumberOfControlPoints,0.0) ;
		for(int i = 0; i<NumberOfControlPoints; i++) {
			PBR.GetBoundaryRepresentation(P,ArcParam[i]) ;
			PointsBoundary[2*i  ] = P[0] ;
			PointsBoundary[2*i+1] = P[1] ;
		}
		return PointsBoundary ;
	}
} ;

#endif