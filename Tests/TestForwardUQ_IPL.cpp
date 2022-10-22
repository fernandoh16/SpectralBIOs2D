#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <complex>

// Sources
#include "Point/Point.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"
#include "BIOsProblems/DirichletProblem.hpp"

// Boundary Representation
#include "BoundaryRepresentations/Boomerang.hpp"

//#include <boost/mpi.hpp>
/*
#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

//using namespace boost ;
using namespace boost::program_options ;
using namespace std ;
*/

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

//using namespace boost ;
using namespace boost::program_options ;
using namespace std ;

// JSON output
#include <jsoncons/json.hpp>

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		return (PBR->NormOfTangentVector(t)) * E * exp(complex<double>(0.0,WaveNumber * (Direction & P))) ;
	}
} ;

complex<double> GetPlaneWave(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->Wave(t) ;
}

// === Quantity of Interest ===
struct QuantityOfInterest {
	ParametricBoundaryRepresentation & PBR ;
	int      NumberOfControlPoints ;
	vector <double>       ArcParam ;
	vector < complex < double > > SolutionBoundary ;
	QuantityOfInterest(ParametricBoundaryRepresentation & PBR_, int NumberOfControlPoints_): 
					PBR(PBR_) {
		NumberOfControlPoints = NumberOfControlPoints_ ;
		ArcParam.assign(NumberOfControlPoints,0.0) ;
		SolutionBoundary.assign(NumberOfControlPoints,complex<double>(0.0,0.0)) ;
		for(int i = 0; i<NumberOfControlPoints; i++) {
			ArcParam[i] = (double)(i)/(double)(NumberOfControlPoints) ;
		}
	}
	vector < complex < double > > ComputeQoI(SpectralQuantity & SQ, vector < double > & y) {
		PBR.SetParameters(y) ;		
		for(int i = 0; i<NumberOfControlPoints; i++) {
			SolutionBoundary[i] = SQ.GetScaled(ArcParam[i],PBR.NormOfTangentVector(ArcParam[i])) ;
		}
		return SolutionBoundary ;
	}
} ;

// === Forward UQ SL BIO ===
int main(int argc, char* argv[]) {
	// === Communicator ===
	/*
    boost::mpi::environment env(argv);
    boost::mpi::communicator world;
    */
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
	int Level ;
	//
	string DGVP = "/Users/FernandoH/Documents/ETH/Codes/UQ_BIOs/SPODWeights" ;
	//
	bool verbose ;  
    options_description desc("SL and HOQMC for the Single Layer BIO in a C^2-boundary.") ;
    desc.add_options()("help,h", "Show this help message")
        ("Sampler"            , value<string>(&Sampler)->default_value("Halton")          , "Sampler (Halton, IPL)")
        ("NumberOfParameters" , value<int>(&NumberOfParameters)->default_value(10)        , "Number Of Parameters (integration Dimensions)") 
        ("NumberOfModes"      , value<int>(&NumberOfModes)->default_value(50)             , "Number Of Spectral Modes") 
        ("NumberOfPoints"     , value<int>(&NumberOfPoints)->default_value(24)            , "Number of Quad points for PBR") 
        ("NumberOfCycles"     , value<int>(&NumberOfCycles)->default_value(1)             , "Number of Quad cycles for PBR") 
        ("GenVectorPath"      , value<string>(&GenVectorPath)->default_value(DGVP)        , "Path to generating vectirs") 
        ("Problem"            , value<string>(&Problem)->default_value("Helmholtz")       , "Problem (Laplace,Helmholtz)") 
        ("WaveNumber"         , value<double>(&WaveNumber)->default_value(1.0)            , "Wave Number")
        ("Magnitude"          , value<double>(&Magnitude)->default_value(1.0)             , "Magnitude of the Incident Field") 
        ("nSamples"           , value<int>   (&nSamples)->default_value(99)               , "Number Of FFT Samples") 
        ("DirectionX"         , value<double>(&DirectionX)->default_value(1.0)            , "X-component of the incident field") 
        ("DirectionY"         , value<double>(&DirectionY)->default_value(0.0)            , "Y-component of the incident field")
        ("NControlPoints"     , value<int>    (&NumberOfControlPoints)->default_value(50) , "Number of Boundary Points") 
        ("Level"              , value<int>    (&Level)->default_value(0)                  , "Level") 
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
	ParametricBoundaryRepresentation PBR(NumberOfParameters,Boomerang,Boomerang_der) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.0) ; 
	PBR.SetParameters(y) ;
	// Quadrature
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	// === Spectral BIOs ===
	SpectralBIOs S(NumberOfModes,Problem,WaveNumber,&PBR) ;
	S.BuildSpectralBIOs(nSamples) ; 
	// === Right Hand Side === 
	SpectralQuantity RHS(NumberOfModes,nSamples) ;
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
	// === Dirichlet Problem ===
	DirichletProblem DP(S,RHS) ;
	DP.BuildMatrix() ;
	DP.BuildRHSIndirectMethod() ;
	// === Get Solution ===
	SpectralQuantity Solution(NumberOfModes) ;
	// DP.GetSolution(Solution) ;
	// === QoI ===
	QuantityOfInterest QoI(PBR,NumberOfControlPoints) ;
    // === Integrand for forward UQ === 
    auto Integrand = [&DP,&Solution,&PBR,&QoI,&PW,&RHS](std::vector<double> y, int level){
    	PBR.SetParameters(y) ;
    	RHS.ConvertToSpectral(GetPlaneWave,&PW) ;
    	DP.BuildMatrix() ;
    	DP.BuildRHSIndirectMethod() ;
    	DP.Solve() ;
    	DP.GetSolution(Solution) ;
    	auto Obs = QoI.ComputeQoI(Solution,y) ;
    	return Obs ;
    } ;
    /*
    auto Obs = Integrand(y,0) ;
	double EvalPoint ;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		EvalPoint = (double)(i)/(double)(NumberOfControlPoints) ;
		cout << (Obs[i]).real() << " " ;
		cout << (Obs[i]).imag() << "\n" ;
	}
	*/
    // === Samples ===
    auto LeveltoFile = [&GenVectorPath](int l) {
        stringstream ss;
        int m = l + 1 ;  // number of points: N=2^m
        ss << GenVectorPath << "/standard_spod_a3_C0.1_SPOD_TZ_t0.2_z3/SPOD_t2.000e-01_z3.000e+00_m" << m << ".json" ;
        return ss.str();
    } ;

    // Starting from Level = 0,1,...
    gMLQMC::samplers::IPL<double> IPL_Sampler(LeveltoFile,NumberOfParameters) ;
    int M = pow(2, Level + 1);
    cout <<  M << "\n" ; 
    auto Out = gMLQMC::SLQMC<gMLQMC::strategy::serial>(Integrand,IPL_Sampler, M, Level) ;
    for(int i = 0; i<NumberOfControlPoints; i++) {
		cout << (Out[i]).real() << "    " << (Out[i]).imag() << "\n" ;
	}
	return 0 ;
}

