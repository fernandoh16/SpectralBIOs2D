#ifndef SpectralBIOs_HEADERDEF 
#define SpectralBIOs_HEADERDEF

#include <cmath> 
#include <complex>
#include <map>
#include <vector>
#include <string>

#include <fftw3.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std ;

#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

class SpectralBIOs {

private:
	//
	int nModes ;
	ParametricBoundaryRepresentation * PBR ;
	int Samples ; 
	int NewSamples ;
	//
	string CaseWaveNumber ;
	double WaveNumber ;
	map < pair< int , int > , complex < double > > V  ;
	map < pair< int , int > , complex < double > > W  ;
	map < pair< int , int > , complex < double > > K  ;
	map < pair< int , int > , complex < double > > Kt ;
	map < pair< int , int > , complex < double > > M ;
	// FFT
	fftw_plan PlanSpectralBIOs ;
	fftw_plan PlanSpectralBIOs_1D ;
	fftw_complex  *  Input     ;
  	fftw_complex  * Output     ;
	fftw_complex  *  Input_1D  ;
  	fftw_complex  * Output_1D  ;
  	// Cross Interactions
	ParametricBoundaryRepresentation * PBR_C ;
	map < int , map < pair< int , int > , complex < double > > > T    ;
	// Is built?
	

public:
	// Constructor
	SpectralBIOs(int nModes_, 
				 ParametricBoundaryRepresentation * PBR_, 
			     int Samples_,
				 string CaseWavenumber_,
				 double WaveNumber_) ;
	// Get Parameters
	int GetNumberofModes()     ;
	int GetNumberOfSamples()   ;
	ParametricBoundaryRepresentation * GetPBR() ;
	// Build BIOs
	void BuildSpectralBIOs()   ;
	void SampleBIO(complex < double > (*f1)(double t, double s, void * Input), 
				   complex < double > (*f2)(double t, void * Input), 
				   void * Input_) ;
	void AssembleBIO(map < pair< int , int > , complex < double > > * BIO_M) ;
	void BuildSpectralV()      ;
	void BuildSpectralK()      ;
	void BuildSpectralV2()     ;
	void BuildSpectralK2()     ;
	void BuildSpectralKt()     ;
	void BuildSpectralW()      ;
	void BuildSpectralW2()     ;
	// Mass Matrix
	void BuildMassMatrix()     ;
	// Cross Interactions
	void SetBoundaryCrossInteraction(ParametricBoundaryRepresentation * PBR_) ;
	void BuildCrossInteraction() ;
	// Get Matrices
	map < pair< int , int > , complex < double > > * GetV()  ; 
	map < pair< int , int > , complex < double > > * GetW()  ; 
	map < pair< int , int > , complex < double > > * GetK()  ;
	map < pair< int , int > , complex < double > > * GetKt() ; 
	map < pair< int , int > , complex < double > > * GetT(int B) ;
	// Destructor
	~SpectralBIOs() ;
};
	
#endif
