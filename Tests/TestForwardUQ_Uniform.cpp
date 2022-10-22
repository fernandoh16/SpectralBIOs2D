#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std ;

// Sources
#include "Point/Point.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "BIOsProblems/DirichletProblem.hpp"

// Boundary Representation
#include "BoundaryRepresentations/Boomerang.hpp"

// gMLQMC
#define GMLQMC_SERIAL_FLAG  // don't include code using MPI or boost::mpi. This needs to be before gMLQMC includes.
#include <gMLQMC/gMLQMC.hpp>
#include <gMLQMC/samplers/uniform.hpp>
#include <gMLQMC/samplers/IPL.hpp>
#include <gMLQMC/strategy/serial.hpp>

// JSON output
// #include <jsoncons/json.hpp>

// Boost
#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options.hpp>

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		return E * exp(complex<double>(0.0,WaveNumber * (Direction & P))) ;
	}
} ;

complex<double> GetPlaneWave(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->Wave(t) ;
}

/*
SpectralQuantity operator*(double a, SpectralQuantity const & S) {
	return S*a ;
}
*/

// === Forward UQ SL BIO ===

int main(int argc, char* argv[]) {
	// === Parametric Boundary Representation ===
	int NumberOfParameters = 10 ;
	int nModes = 10 ;
	ParametricBoundaryRepresentation PBR(NumberOfParameters,Boomerang,Boomerang_der) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.0) ; 
	PBR.SetParameters(y) ;
	// Quadrature
	int NumberOfPoints = 24 ;
	int NumberOfCycles =  1 ; 
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	// === Spectral BIOs ===
	double WaveNumber = 1.0 ;
	double E = 1.0 ;
	int nSamples = 99 ;
	SpectralBIOs S(nModes,"Helmholtz",WaveNumber,&PBR) ;
	S.BuildSpectralBIOs(nSamples) ; 
	// === Right Hand Side === 
	SpectralQuantity RHS(nModes,nSamples) ;
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = E ;
	PW.Direction = Point(1.0,0.0) ;
	PW.PBR = &PBR ;
	RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
	// === Dirichlet Problem ===
	DirichletProblem DP(S,RHS) ;
	DP.BuildMatrix() ;
	DP.BuildRHSIndirectMethod() ;
	// === Get Solution ===
	SpectralQuantity Solution(nModes) ;
	DP.GetSolution(Solution) ;
    // === Integrand for forward UQ === 
    auto Integrand = [&DP,&Solution](std::vector<double> y, int level){
    	DP.Solve() ;
    	DP.GetSolution(Solution) ;
        return Solution ;
    } ;
    // === Samples ===
    gMLQMC::samplers::Uniform<double> Sampler(NumberOfParameters);

    // === Approximate the Expectation ===
    int M = 1.0 ;
    int L = 1.0 ;
    auto Out = gMLQMC::SLQMC<gMLQMC::strategy::serial>(Integrand,Sampler, M, L) ;
    // === Print Solution ===
	int NumberOfPointsPrint = 100 ;
	double EvalPoint ;
	for(int i = 0; i<NumberOfPointsPrint; i++) {
		EvalPoint = (double)(i)/(double)(NumberOfPointsPrint) ;
		cout << Out.GetScaled(EvalPoint,PBR.NormOfTangentVector(EvalPoint)) << "\n" ;
	}
	return 0 ;
}

