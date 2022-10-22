#ifndef PLANE_WAVE_HPP
#define PLANE_WAVE_HPP
#pragma once

using namespace std ;

#include "Point/Point.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		// return (PBR->NormOfTangentVector(t)) * E * exp(complex<double>(0.0,WaveNumber * (Direction & P))) ;
		return -E * exp(complex<double>(0.0,WaveNumber*(Direction & P))) ;
	}
} ;

complex<double> GetPlaneWave(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->Wave(t) ;
}

complex<double> GetPlaneWaveNormalized(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	double Norm = (PW->PBR)->NormOfTangentVector(t) ;
	return (PW->Wave(t)) * Norm ;
}

#endif