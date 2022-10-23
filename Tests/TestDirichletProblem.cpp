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
/*
#include "BoundaryRepresentations/Boomerang.hpp"
#include "BoundaryRepresentations/Kidney.hpp"
#include "BoundaryRepresentations/Circle.hpp"
#include "BoundaryRepresentations/Fourier.hpp"
#include "BoundaryRepresentations/CCoeffAk.hpp" 
*/
#include "BoundaryRepresentations/Circle.hpp"

#include <boost/algorithm/string.hpp>  // to_upper
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "Halton/halton.hpp"
#include "Halton/halton.cpp"

// #include "BayesianInversion/InverseScatteringBI/PlaneWave.hpp"

struct PlaneWave {
	double WaveNumber ;
	double E ;
	Point Direction ;
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	Point N ;
	Point P_ext;
	complex< double > GetDirichletTrace(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		//return -1.0/(2*M_PI)*log(length(P,P_ext));
		return exp(complex<double>(0.0,WaveNumber*(Direction & P))) ;
	} 
	complex< double > GetNeumannTrace(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		PBR->GetNormalVector(N,t);
		// return -1.0/(2*M_PI)*(N&(P-P_ext))/pow(length(P,P_ext),2.0);
		return complex < double > (0.0,WaveNumber*(N & Direction))*exp(complex<double>(0.0,WaveNumber*(Direction & P))) ;
	} 
	complex< double > GetNeumannTraceScaled(double t) {	
		PBR->GetBoundaryRepresentation(P,t) ;
		PBR->GetNormalVector(N,t);
		// return -1.0/(2*M_PI)*(N&(P-P_ext))/pow(length(P,P_ext),2.0);
		return complex < double > (0.0,WaveNumber*(N & Direction))*exp(complex<double>(0.0,WaveNumber*(Direction & P)))*(PBR->NormOfTangentVector(t)) ;
	} 
} ;

complex < double > GetDirichletTrace(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->GetDirichletTrace(t) ;
}

complex < double > GetNeumannTrace(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->GetNeumannTrace(t) ;
}

complex < double > GetNeumannTraceScaled(double t, void * Input) {
	PlaneWave * PW = (PlaneWave *) Input ;
	return PW->GetNeumannTraceScaled(t) ;
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
	double ZetaDecay ;
	double Sigma ;
	//
	bool verbose ;  
    options_description desc("-- Test Dirichlet Problem --") ;
    desc.add_options()("help,h", "Show this help message")
        ("NumberOfParameters" , value<int>(&NumberOfParameters)->default_value(1)       , "Number Of Parameters (integration Dimensions)") 
        ("NumberOfModes"      , value<int>(&NumberOfModes)->default_value(50)             , "Number Of Spectral Modes") 
        ("NumberOfPoints"     , value<int>(&NumberOfPoints)->default_value(24)            , "Number of Quad points for PBR") 
        ("NumberOfCycles"     , value<int>(&NumberOfCycles)->default_value(10)            , "Number of Quad cycles for PBR") 
        ("Problem"            , value<string>(&Problem)->default_value("Laplace")       , "Problem (Laplace,Helmholtz)") 
        ("WaveNumber"         , value<double>(&WaveNumber)->default_value(0.0)            , "Wave Number")
        ("Magnitude"          , value<double>(&Magnitude)->default_value(1.0)             , "Magnitude of the Incident Field") 
        ("nSamples"           , value<int>   (&nSamples)->default_value(150)              , "Number Of FFT Samples") 
        ("DirectionX"         , value<double>(&DirectionX)->default_value(1.0)            , "X-component of the incident field") 
        ("DirectionY"         , value<double>(&DirectionY)->default_value(0.0)            , "Y-component of the incident field")
        ("NControlPoints"     , value<int>    (&NumberOfControlPoints)->default_value(50) , "Number of Boundary Points") 
        ("Solution"           , value<string>(&PrintSolution)->default_value("solution.txt")   , "Solution") 
		("ZetaDecay"          , value<double>(&ZetaDecay)->default_value(2.0)                       , "Decay Fourier Series")
		("Sigma"              , value<double>(&Sigma)->default_value(0.15)                          , "Sigma")
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
	// Parametric Boundary Representation
	/*
	double Index ;
	vector < double > Coeff_A ;
	for(int i = 1; i<=NumberOfParameters; i++) {
		if(i % 2 == 0) {
			Index = (double)(i)/2.0 ;
		}
		else {
			Index = (double)(i+1)/2.0 ;
		}
		Coeff_A.push_back(pow(Index,-ZetaDecay)) ;
	}
	ParametricBoundaryRepresentation PBR(NumberOfParameters,Coeff_A,Sigma) ;
	PBR.SetParametricBoundaryRepresentation(CircleRep2,CircleRep2Der,CircleRep2SecDer) ;
	PBR.SetParametricBoundaryRepresentationBasis(FourierPBR,FourierPBRDer,FourierPBRSecDer) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.5) ; 
	PBR.SetParameters(y) ;
	*/
	// Parametric Boundary Representation
	ParametricBoundaryRepresentation PBR(NumberOfParameters,R_Fourier,R_der_Fourier,R_sec_der_Fourier) ;
	vector <double> y ;
	y.assign(NumberOfParameters,0.5) ; 
	PBR.SetParameters(y) ;
	// Quadrature 
	PBR.SetQuadrature(NumberOfPoints,NumberOfCycles) ;
	// Spectral BIOs 
	SpectralBIOs S(NumberOfModes,&PBR,nSamples,Problem,WaveNumber) ;
	// S.BuildSpectralBIOs() ;
	// S.BuildSpectralV();
	// S.BuildSpectralV2();
	//
	S.BuildSpectralV2();
	S.BuildSpectralK2();
	//  Right Hand Side  
	SpectralQuantity RHS1(NumberOfModes,nSamples) ;
	SpectralQuantity RHS2(NumberOfModes,nSamples) ;
	// Plane Wave
	PlaneWave PW ;
	PW.WaveNumber = WaveNumber ;
	PW.E = Magnitude ;
	PW.Direction = Point(DirectionX,DirectionY) ;
	PW.PBR = &PBR ;
	PW.P_ext = Point(1.5,0.0,0.0);
	//
	RHS1.ConvertToSpectral(GetDirichletTrace,&PW) ;
	RHS2.ConvertToSpectral(GetNeumannTraceScaled,&PW) ;
	// Dirichlet Problem
	DirichletProblem DP(S) ;
	DP.BuildMatrix() ;
	DP.BuildRHSDirectMethod(RHS1) ;
	DP.Solve() ;
	// Get Solution
	SpectralQuantity Solution(NumberOfModes) ;
	DP.GetSolution(Solution) ;
	// Print Solution
	ofstream print_output(PrintSolution) ;
	assert(print_output.is_open()) ;
	// Post-Processing
	double EvalPoint ;
	complex < double > Value_Exact;
	for(int i = 0; i<NumberOfControlPoints; i++) {
		EvalPoint = (double)(i)/(double)(NumberOfControlPoints) ;
		print_output << EvalPoint << " " ;
		Value_Exact = GetNeumannTrace(EvalPoint,&PW);
		print_output << setprecision(10) << (Solution.GetScaled(EvalPoint,PBR.NormOfTangentVector(EvalPoint))).real()<< " " << Value_Exact.real() << "\n" ;
	}
	print_output.close();
	double Error = 0;
	for(int i = -NumberOfModes; i<=NumberOfModes; i++) {
		Error = Error + norm(Solution[i]-RHS2[i]);
	}
	/*
	cout << "Error:" << sqrt(Error) << "\n";
	// Halton Points
	int n_points = 10;
	int dim = 10;
	double * I;
	I = halton(10,dim);
	for (int i=0; i<dim; i++) {
		cout << I[i] << "\n";
	}
	cout << "Hola" << "\n";
	*/
	return 0 ;
}

