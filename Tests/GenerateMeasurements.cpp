#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>
#include <random>

using namespace std ;

// Sources
#include "Point/Point.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "BIOsProblems/DirichletProblem.hpp"
#include "BoundaryPotentials/BoundaryPotentials.hpp"

// Bayesian Inversion
#include "BayesianInversion/BayesianInversion.hpp"
#include "BayesianInversion/InverseScatteringBI/ForwardModel.hpp"
#include "BayesianInversion/InverseScatteringBI/PlaneWave.hpp"
#include "BayesianInversion/InverseScatteringBI/Potential.hpp"
#include "BayesianInversion/InverseScatteringBI/QuantityOfInterest.hpp"

// Boundary Representation
#include "BoundaryRepresentations/Kidney.hpp"

#include <direct_product.hpp>
#include <gMLQMC/tools/vector_ops.hpp>

#define GMLQMC_SERIAL_FLAG  // don't include code using MPI or boost::mpi. This needs to be before gMLQMC includes.
#include <gMLQMC/gMLQMC.hpp>
#include <gMLQMC/samplers/uniform.hpp>
#include <gMLQMC/samplers/IPL.hpp>
#include <gMLQMC/strategy/serial.hpp>
#include <gMLQMC/samplers/halton.hpp>

#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace boost::program_options ;
using namespace std ;

// JSON I/O
#include <jsoncons/json.hpp>

// ==================================================
// === Bayesian Inversion for Acoustic Scattering ===
// ==================================================

int main(int argc, char* argv[]) {
	// === Inputs ===
	string Sampler ;
	int NumberOfParameters    ;
	int NumberOfModes         ;
	int NumberOfPoints        ;
	int NumberOfCycles        ;
	string GenVectorPath      ;
	string Problem            ;
	double WaveNumber         ;
	double Magnitude          ;
	int nSamples              ;
	double DirectionX         ;
	double DirectionY         ;
	int NumberOfControlPoints ;
	int NumberOfDirections    ;
	int Level ;
	double Gamma ;
	//
	string DGVP = "/Users/FernandoH/Documents/ETH/Codes/UQ_BIOs/SPODWeights" ;
	//
	string MeasFile ;
	string TrueBoundary ;
	//
	bool verbose ;  
    options_description desc("--GM--") ;
    desc.add_options()("help,h", "Show this help message")
        ("Sampler"            , value<string>(&Sampler)->default_value("Halton")          , "Sampler (Halton, IPL)")
        ("NumberOfParameters" , value<int>(&NumberOfParameters)->default_value(100)       , "Number Of Parameters (integration Dimensions)") 
        ("NumberOfModes"      , value<int>(&NumberOfModes)->default_value(25)             , "Number Of Spectral Modes") 
        ("NumberOfPoints"     , value<int>(&NumberOfPoints)->default_value(24)            , "Number of Quad points for PBR") 
        ("NumberOfCycles"     , value<int>(&NumberOfCycles)->default_value(2)             , "Number of Quad cycles for PBR") 
        ("GenVectorPath"      , value<string>(&GenVectorPath)->default_value(DGVP)        , "Path to generating vectors") 
        ("Problem"            , value<string>(&Problem)->default_value("Helmholtz")       , "Problem (Laplace,Helmholtz)") 
        ("WaveNumber"         , value<double>(&WaveNumber)->default_value(1.0)            , "Wave Number")
        ("Magnitude"          , value<double>(&Magnitude)->default_value(1.0)             , "Magnitude of the Incident Field") 
        ("Gamma"              , value<double>(&Gamma)->default_value(1.0)                 , "Gamma") 
        ("nSamples"           , value<int>   (&nSamples)->default_value(50)              , "Number Of FFT Samples") 
        ("DirectionX"         , value<double>(&DirectionX)->default_value(1.0)            , "X-component of the incident field") 
        ("DirectionY"         , value<double>(&DirectionY)->default_value(0.0)            , "Y-component of the incident field")
        ("NControlPoints"     , value<int>    (&NumberOfControlPoints)->default_value(50) , "Number of Boundary Points") 
        ("Level"              , value<int>    (&Level)->default_value(0)                  , "Level") 
        ("NDirections"        , value<int>(&NumberOfDirections)->default_value(16)        , "Number of Observations Directions") 
        ("MeasFile"           , value<string>(&MeasFile)->default_value("MeasFile.txt")   , "Measurements Output File") 
        ("TrueBoundary"       , value<string>(&TrueBoundary)->default_value("TrueBoundary.txt")   , "True Boundary Output File") 
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
	ParametricBoundaryRepresentation PBR(NumberOfParameters,Kidney_rep,Kidney_rep_der) ;
	vector <double> y ;
	//y.assign(NumberOfParameters,0.8127) ; 
	y.assign(NumberOfParameters,0.5) ; 
	PBR.SetParameters(y) ;
	// Quadrature
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	
	// === Spectral BIOs ===
	SpectralBIOs S(NumberOfModes,Problem,WaveNumber,&PBR,nSamples) ;
	S.BuildSpectralV() ; 
	// === Right Hand Side === 
	SpectralQuantity RHS(NumberOfModes,nSamples) ;
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
	
	// === Dirichlet Problem ===
	DirichletProblem DP(S) ;
	DP.BuildMatrix() ;
	DP.BuildRHSIndirectMethod(RHS) ;
	
	// === Get Solution ===
	SpectralQuantity Solution(NumberOfModes) ;
	DP.GetSolution(Solution) ;
	
	// === Boundary Potentials & Far Field ===
	BoundaryPotentials BP(Solution,PBR,nSamples,WaveNumber) ;
	ForwardModel FM(PBR,DP,BP,RHS,Solution,PW,S,NumberOfDirections) ;

	// === Generate Measurements ===
	vector < double > y_star ;
	y_star.assign(NumberOfParameters,1.0) ;
	vector < complex < double > > M ;
	M.assign(NumberOfDirections * NumberOfDirections, complex < double > (0.0,0.0)) ; 
	M = FM(y_star) ;
	// Perturbation
	random_device rd ;
	std::mt19937 gen(rd());  
	normal_distribution<double> distribution(0.0,Gamma) ;
	for(int i = 0; i<M.size(); i++) {
		M[i] = M[i] + complex<double>(distribution(gen),distribution(gen)) ;
	}
	// Print Measurements
	ofstream Output ;
	Output.open(MeasFile) ;
	for(int i = 0; i<M.size(); i++) {
		Output << M[i].real() << "   " << M[i].imag()<< "\n" ;
	}
	Output.close() ;
	
	// === QoI ===
	QuantityOfInterest QoI(PBR,NumberOfControlPoints) ;
	
    // === Print Output ==
    ofstream OutputBoundary ;
    //
    OutputBoundary.open(TrueBoundary) ;
    //
    auto Out = QoI(y_star) ;
    for(int i = 0; i<NumberOfControlPoints; i++) {
    	OutputBoundary << Out[2 * i] << "   " << Out[2 * i + 1]<< "\n" ; 
   	}
    OutputBoundary.close();
	return 0 ;
}

