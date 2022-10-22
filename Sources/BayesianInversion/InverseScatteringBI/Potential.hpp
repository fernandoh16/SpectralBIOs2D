#ifndef POTENTIAL_HPP
#define POTENTIAL_HPP
#pragma once

using namespace std ;

#include "Point/Point.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "BIOsProblems/DirichletProblem.hpp"
#include "BoundaryPotentials/BoundaryPotentials.hpp"

struct Potential {
	double Gamma ;
	double operator()(const vector < complex < double > > & M1, const vector < complex < double > >  & M2) const {
		double Out  = 0;
		assert(M1.size() == M2.size()) ;
		for(int i = 0; i<M1.size(); i++) {
			Out = Out + pow(abs(M1[i]-M2[i]),2.0) ;
		} 
		return 0.5 * Out / Gamma ;
	} 
} ;

#endif