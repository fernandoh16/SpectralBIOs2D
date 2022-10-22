#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
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
// #include "BoundaryRepresentations/Boomerang.hpp"
#include "BoundaryRepresentations/Circle.hpp"

#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

// #include "BayesianInversion/InverseScatteringBI/PlaneWave.hpp"

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		return 1+0.5*sin(2.0*M_PI*t) + pow(0.5,4)*sin(3.0*2.0*M_PI*t) + pow(0.5,4)*cos(4.0*2.0*M_PI*t) ;
		// return  -E * exp(complex<double>(0.0,WaveNumber*(Direction & P))) ;
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

using namespace boost::program_options ;
using namespace std ;

int main(int argc, char* argv[]) {
	// === Inputs ===
	int NumberOfParameters    ;
	int NumberOfModes         ;
	int NumberOfPoints        ;
	int NumberOfCycles        ;
	string Problem            ;
	double WaveNumber         ;
	double Magnitude          ;
	int nSamples              ;
	double DirectionX         ;
	double DirectionY         ;
	int NumberOfControlPoints ;
	//
	string PrintSolution ;
	//
	bool verbose ;  
    options_description desc("-- Test Dirichlet Problem --") ;
    desc.add_options()("help,h", "Show this help message")
        ("NumberOfParameters" , value<int>(&NumberOfParameters)->default_value(1)       , "Number Of Parameters (integration Dimensions)") 
        ("NumberOfModes"      , value<int>(&NumberOfModes)->default_value(50)             , "Number Of Spectral Modes") 
        ("NumberOfPoints"     , value<int>(&NumberOfPoints)->default_value(24)            , "Number of Quad points for PBR") 
        ("NumberOfCycles"     , value<int>(&NumberOfCycles)->default_value(10)            , "Number of Quad cycles for PBR") 
        ("Problem"            , value<string>(&Problem)->default_value("Helmholtz")       , "Problem (Laplace,Helmholtz)") 
        ("WaveNumber"         , value<double>(&WaveNumber)->default_value(1.0)            , "Wave Number")
        ("Magnitude"          , value<double>(&Magnitude)->default_value(1.0)             , "Magnitude of the Incident Field") 
        ("nSamples"           , value<int>   (&nSamples)->default_value(350)              , "Number Of FFT Samples") 
        ("DirectionX"         , value<double>(&DirectionX)->default_value(1.0)            , "X-component of the incident field") 
        ("DirectionY"         , value<double>(&DirectionY)->default_value(0.0)            , "Y-component of the incident field")
        ("NControlPoints"     , value<int>    (&NumberOfControlPoints)->default_value(50) , "Number of Boundary Points") 
        ("Solution"           , value<string>(&PrintSolution)->default_value("solution.txt")   , "Solution") 
        ("verbose,v"          , value<bool>(&verbose)->default_value(false)->implicit_value(true)   , "Verbose output")
    ;
    variables_map vm;
    try {
        store(parse_command_line(argc, argv, desc), vm);
    } catch(...) {
        cout << desc << std::endl;
        return 1;
    }
    notify(vm);
    // print help if --help or -h specified
    if (vm.count("help") || vm.count("h")) {
        cout << desc << std::endl;
        return 1;
    } 
	// === Parametric Boundary Representation ===
	ParametricBoundaryRepresentation PBR(NumberOfParameters,R_Fourier,R_der_Fourier,R_sec_der_Fourier) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.5) ; 
	PBR.SetParameters(y) ;
	// Quadrature 
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	// === Spectral BIOs ===
	SpectralBIOs S(NumberOfModes,&PBR,nSamples,Problem,WaveNumber) ;
	S.BuildSpectralBIOs() ; 
	// === Right Hand Side === 
	SpectralQuantity RHS1(NumberOfModes,nSamples) ;
	SpectralQuantity RHS2(NumberOfModes,nSamples) ;
	//
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	//
	RHS1.ConvertToSpectral(GetPlaneWave,&PW) ;
	RHS2.ConvertToSpectral(GetPlaneWaveNormalized,&PW) ;
	// === Dirichlet Problem ===
	DirichletProblem DP(S) ;
	DP.BuildMatrix() ;
	DP.BuildRHSDirectMethod(RHS1,RHS2) ;
	// DP.BuildRHSIndirectMethod() ;
	DP.Solve() ;
	// === Get Solution ===
	SpectralQuantity Solution(NumberOfModes) ;
	DP.GetSolution(Solution) ;
	// === Print Matrix ===
	ofstream print_output(PrintSolution) ;
	assert(print_output.is_open()) ;
	//
	double EvalPoint ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		EvalPoint = (double)(i)/(double)(NumberOfControlPoints) ;
		print_output << EvalPoint << " " ;
		/*
		//print_output << (Solution(EvalPoint)).real() << " " ;
		//print_output << (Solution(EvalPoint)).imag() << "\n" ;
		*/
		print_output << (Solution.GetScaled(EvalPoint,PBR.NormOfTangentVector(EvalPoint))).real() << "  " ;
		print_output << (Solution.GetScaled(EvalPoint,PBR.NormOfTangentVector(EvalPoint))).imag() << "\n" ;
	}
	return 0 ;
}

