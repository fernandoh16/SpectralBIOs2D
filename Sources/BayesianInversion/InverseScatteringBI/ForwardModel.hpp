#ifndef FORWARD_MODEL_HPP
#define FORWARD_MODEL_HPP
#pragma once

using namespace std ;

#include "Point/Point.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "BIOsProblems/DirichletProblem.hpp"
#include "BoundaryPotentials/BoundaryPotentials.hpp"
#include "BayesianInversion/InverseScatteringBI/PlaneWave.hpp"

struct ForwardModel {
	//
	ParametricBoundaryRepresentation & PBR ;
	DirichletProblem & DP ;
	BoundaryPotentials & BP ;
	SpectralQuantity & RHS ;
	SpectralQuantity & Sol ;
	PlaneWave & PW ;
	SpectralBIOs & S ;
	//
	vector< complex < double > > Ms ;
	vector < Point * >    IncidentDirections ;
	vector < Point * > ObservationDirections ;
	//
	ForwardModel(ParametricBoundaryRepresentation & PBR_, 
				 DirichletProblem & DP_, 
				 BoundaryPotentials & BP_,  
				 SpectralQuantity & RHS_,
				 SpectralQuantity & Sol_,
				 PlaneWave & PW_,
				 SpectralBIOs & S_,
				 int NumberOfDirections) :
		PBR(PBR_)      ,
		DP(DP_)        ,
		BP(BP_)        ,
		RHS(RHS_)      ,
		Sol(Sol_)      ,      
		PW(PW_)        ,
		S(S_)          {
		double x, y ;
		for(int i = 0; i<NumberOfDirections; i++) {
			x = cos(2.0 * M_PI *(double)(i)/(double)(NumberOfDirections)) ;
			y = sin(2.0 * M_PI *(double)(i)/(double)(NumberOfDirections)) ;
			IncidentDirections.push_back(new Point(x,y)) ;
			ObservationDirections.push_back(new Point(x,y)) ;
		}
		Ms.assign(NumberOfDirections * NumberOfDirections, complex<double>(0.0,0.0)) ;
	}
	vector< complex < double > > & operator()(vector < double > & y) {
		for(int i = 0; i<IncidentDirections.size(); i++) {
			PBR.SetParameters(y) ;
			PW.Direction = *(IncidentDirections[i]) ;
			RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
			S.BuildSpectralBIOs() ; 
			DP.BuildMatrix() ;
			DP.BuildRHSIndirectMethod(RHS) ;
   			DP.Solve() ;
    		DP.GetSolution(Sol) ;
    		for(int j = 0; j<ObservationDirections.size() ; j++) {
    			Ms[j+(ObservationDirections.size())*i] = BP.FarFieldSingleLayerPotential(*(ObservationDirections[j])) ;
    		}
		}
		return Ms ; 
	}
} ;

#endif