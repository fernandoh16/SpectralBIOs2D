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

// #define GMLQMC_SERIAL_FLAG  // don't include code using MPI or boost::mpi. This needs to be before gMLQMC includes.
#include <gMLQMC/gMLQMC.hpp>
#include <gMLQMC/samplers/uniform.hpp>
#include <gMLQMC/samplers/IPL.hpp>
#include <gMLQMC/strategy/full.hpp>
#include <gMLQMC/samplers/halton.hpp>

#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>

using namespace boost::program_options ;
using namespace std ;

// JSON I/O
#include <jsoncons/json.hpp>

#define pout if (world.rank() == 0) std::cout

// ==================================================
// === Bayesian Inversion for Acoustic Scattering ===
// ==================================================

int main(int argc, char* argv[]) {
	// === Communicator ===
    boost::mpi::environment env(argv);
    boost::mpi::communicator world;
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
	string Prior  ;
	string Posterior ;
	string Boundary ;
	//
	string MeasurementsFile ;
	//
	bool verbose ;  
	int IPLzeta ;
    options_description desc("--Bayesian Inversion for the SL BIO--") ;
    desc.add_options()("help,h", "Show this help message")
        ("Sampler"            , value<string>(&Sampler)->default_value("Halton")            , "Sampler (Halton, IPL)")
        ("NumberOfParameters" , value<int>(&NumberOfParameters)->default_value(48)         , "Number Of Parameters (integration Dimensions)") 
        ("NumberOfModes"      , value<int>(&NumberOfModes)->default_value(50)               , "Number Of Spectral Modes") 
        ("NumberOfPoints"     , value<int>(&NumberOfPoints)->default_value(24)              , "Number of Quad points for PBR") 
        ("NumberOfCycles"     , value<int>(&NumberOfCycles)->default_value(10)              , "Number of Quad cycles for PBR") 
        ("GenVectorPath"      , value<string>(&GenVectorPath)->default_value(DGVP)          , "Path to generating vectirs") 
        ("Problem"            , value<string>(&Problem)->default_value("Helmholtz")         , "Problem (Laplace,Helmholtz)") 
        ("WaveNumber"         , value<double>(&WaveNumber)->default_value(1.0)              , "Wave Number")
        ("Magnitude"          , value<double>(&Magnitude)->default_value(1.0)               , "Magnitude of the Incident Field") 
        ("Gamma"              , value<double>(&Gamma)->default_value(1.0)                   , "Magnitude of the Incident Field") 
        ("nSamples"           , value<int>   (&nSamples)->default_value(99)                 , "Number Of FFT Samples") 
        ("DirectionX"         , value<double>(&DirectionX)->default_value(1.0)              , "X-component of the incident field") 
        ("DirectionY"         , value<double>(&DirectionY)->default_value(0.0)              , "Y-component of the incident field")
        ("NControlPoints"     , value<int>    (&NumberOfControlPoints)->default_value(50)   , "Number of Boundary Points") 
        ("Level"              , value<int>    (&Level)->default_value(0)                    , "Level") 
        ("NDirections"        , value<int>(&NumberOfDirections)->default_value(16)          , "Number of Observations Directions") 
        ("Prior"              , value<string>(&Prior)->default_value("prior.txt")           , "Prior") 
        ("Posterior"          , value<string>(&Posterior)->default_value("posterior.txt")   , "Posterior")
        ("Boundary"           , value<string>(&Boundary)->default_value("boundary.txt")     , "Boundary") 
        ("Measurements"       , value<string>(&MeasurementsFile)->default_value("Meas.txt") , "Measurements") 
        ("IPLZeta"            , value<int>(&IPLzeta)->default_value(2)                      , "Number of Quad cycles for PBR") 
        ("verbose,v"          , value<bool>(&verbose)->default_value(false)->implicit_value(true)   , "Verbose output")
    ;
    variables_map vm;
    try {
        store(parse_command_line(argc, argv, desc), vm);
    } catch(...) {
        pout << desc << std::endl;
        return 1;
    }
    notify(vm);
    // print help if --help or -h specified
    if (vm.count("help") || vm.count("h")) {
        pout << desc << std::endl;
        return 1;
    } 
	// === Parametric Boundary Representation ===
	ParametricBoundaryRepresentation PBR(NumberOfParameters,Kidney_rep,Kidney_rep_der) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.5) ; 
	PBR.SetParameters(y) ;
	// Quadrature
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	
	// === Spectral BIOs ===
	SpectralBIOs S(NumberOfModes,Problem,WaveNumber,&PBR,nSamples) ;
	S.BuildSpectralBIOs() ; 
	
	// === Right Hand Side === 
	SpectralQuantity RHS(NumberOfModes,nSamples) ;
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	RHS.ConvertToSpectral(GetPlaneWaveNormalized,&PW) ;
	
	// === Dirichlet Problem ===
	DirichletProblem DP(S) ;
	DP.BuildMatrix() ;
	DP.BuildRHSIndirectMethod(RHS) ;
	
	// === Get Solution ===
	SpectralQuantity Solution(NumberOfModes) ;
	
	// === Boundary Potentials & Far Field ===
	BoundaryPotentials BP(Solution,PBR,nSamples,WaveNumber) ;
	ForwardModel FM(PBR,DP,BP,RHS,Solution,PW,S,NumberOfDirections) ;
	
	// === Measurements ===
	vector < complex < double > > M ;
	M.assign(NumberOfDirections * NumberOfDirections, complex < double > (0.0,0.0)) ; 
	if (world.rank() == 0) {
		ifstream MeasInput(MeasurementsFile) ;
		assert(MeasInput.is_open());
		double Aux[2] ;
		for(int i = 0; i<NumberOfDirections; i++) {
			for(int j =0; j<NumberOfDirections; j++) {
				MeasInput >> Aux[0] >> Aux[1] ;
				M[j+NumberOfDirections*i] = complex<double> (Aux[0],Aux[1]) ;
			}
		}
		MeasInput.close();
	}
	broadcast(world,M,0) ;
	
	// === Potential ===
	Potential P ;
	P.Gamma = Gamma ;
	
	// === QoI ===
	QuantityOfInterest QoI(PBR,NumberOfControlPoints) ;
	
	// === Bayesian Inversion Environment ===
	typedef vector < complex < double > > VC;
	typedef vector < double> VD;
   	BayesianInversion<ForwardModel,Potential,VC,QuantityOfInterest,VD,VD> BI(FM,P,M,QoI) ; 
    // === Integrand for BI === 
    auto Integrand = [&BI,&world](std::vector<double> y, int level){
        return BI(y) ;
    } ;
    // === Samples ===
    auto LeveltoFile = [&GenVectorPath,&IPLzeta](int l) {
        stringstream ss;
        int m = l + 1 ;  // number of points: N=2^m
        // ss << GenVectorPath << "/standard_spod_a4_C0.1_SPOD_TZ_t0.2_z4/SPOD_t2.000e-01_z4.000e+00_m" << m << ".json" ;
        ss << GenVectorPath << "/standard_spod_a" << IPLzeta <<"_C0.1_SPOD_TZ_t0.2_z"<< IPLzeta << "/SPOD_t2.000e-01_z" << IPLzeta <<".000e+00_m" << m << ".json" ;
        return ss.str();
    } ;

    // === HoQMC Integration ===
    long Me = pow(2, Level + 1);
	gMLQMC::samplers::IPL<double> IPL_Sampler(LeveltoFile,NumberOfParameters) ;
	auto Out = gMLQMC::SLQMC<gMLQMC::strategy::full>(Integrand,IPL_Sampler,Me,Level,world.size()) ;

    // === Print Output ===
    ofstream PriorOutput ;
    ofstream PosteriorOutput ;
    //
    PriorOutput.open(Prior) ;
    PosteriorOutput.open(Posterior) ;
    //
    if(world.rank() == 0) {
		VD Exp1 =  get<0>(Out);
    	for(int i = 0; i<NumberOfControlPoints; i++) {
    		PriorOutput << fixed << setprecision(15) << Exp1[2 * i] << "   " << Exp1[2 * i + 1] << "\n" ; 
		}
    	VD Exp2 =  get<1>(Out);
   		for(int i = 0; i<NumberOfControlPoints; i++) {
    		PosteriorOutput << fixed << setprecision(15) << Exp2[2 * i] / get<2>(Out) << "   " << Exp2[2 * i + 1] / get<2>(Out) << "\n" ; 
    	}
    }
    PriorOutput.close();
    PosteriorOutput.close();
	return 0 ;
}

