#ifndef BOUNDARYPOTENTIALS_HPP
#define BOUNDARYPOTENTIALS_HPP

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

// Sources
#include "Point/Point.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

using namespace std ;

class BoundaryPotentials {

private:
	ParametricBoundaryRepresentation * PBR ;
	SpectralQuantity * SQ ;
	int nSamples ;
	double WaveNumber ; 
	
public:
	BoundaryPotentials(SpectralQuantity & SQ_, ParametricBoundaryRepresentation & PBR_, int nSamples_, double WaveNumber_) ;
	complex <double> SingleLayerPotential(Point & P) ;
	complex <double> DoubleLayerPotential(Point & P) ;
	complex <double> FarFieldSingleLayerPotential(Point & P) ;
	complex <double> FarFieldDoubleLayerPotential(Point & P) ;
	
} ;
#endif

				