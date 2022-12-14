#include <cmath> 
#include <complex>
#include <iostream>
#include <map> 
#include <cassert>

// FFT
#include <fftw3.h>

// Boost
#include <boost/math/special_functions/bessel.hpp>  
#include <boost/math/special_functions/hankel.hpp>  

// Sources
#include "SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"
#include "ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp"

using namespace std;

int GetIndex(int i, int S) {
	if(i >=0) {
		return i ;
	}
	else {
		return 2*S + 1 + i ;
	}
}

SpectralBIOs::SpectralBIOs(int nModes_, 
						   ParametricBoundaryRepresentation * PBR_, 
						   int Samples_,
						   string CaseWaveNumber_,
				           double WaveNumber_) {
	nModes = nModes_ ;
	PBR = PBR_ ;
	Samples = Samples_ ;
	NewSamples = 2 * Samples + 1 ;
	//
  	Input  = new fftw_complex[NewSamples * NewSamples] ;
  	Output = new fftw_complex[NewSamples * NewSamples] ;
  	//
  	Input_1D  = new fftw_complex[NewSamples] ;
  	Output_1D = new fftw_complex[NewSamples] ;
  	//
	PlanSpectralBIOs    = fftw_plan_dft_2d(NewSamples,NewSamples,Input,Output,FFTW_FORWARD,FFTW_ESTIMATE) ;
	PlanSpectralBIOs_1D = fftw_plan_dft_1d(NewSamples,Input_1D,Output_1D,FFTW_FORWARD,FFTW_ESTIMATE) ;
	//
	CaseWaveNumber = CaseWaveNumber_ ;
	WaveNumber     =     WaveNumber_ ;
}

int SpectralBIOs::GetNumberofModes() {
	return nModes ; 
}

int SpectralBIOs::GetNumberOfSamples() {
	return Samples ; 
}

ParametricBoundaryRepresentation * SpectralBIOs::GetPBR() {
	return PBR ;
}

void SpectralBIOs::BuildSpectralBIOs() {
	BuildSpectralV() ;
	BuildSpectralK() ;
	//BuildSpectralKt() ;
	//BuildSpectralW() ;
	//BuildCrossInteraction() ;
}

void SpectralBIOs::SampleBIO(complex < double > (*f1)(double t, double s, void * Input), 
							 complex < double > (*f2)(double t, void * Input), 
							 void * Input_) {
    double I1 ; 
  	double I2 ;
  	//
	complex < double > Value;
	//
  	for(int i = 0; i<NewSamples; i++) {
  		I1 = (double)(i)/(double)(NewSamples) ;
		for(int j = 0; j<NewSamples; j++) {
			I2 = (double)(j)/(double)(NewSamples) ;
			if(i != j) {
				Value = (*f1)(I1,I2,Input_) ;
				Input[j + NewSamples * i][0] = Value.real() ;
  				Input[j + NewSamples * i][1] = Value.imag() ; 
  			}
  			if (i == j) {
				Value = (*f2)(I1,Input_) ;
  				Input[j + NewSamples * i][0] = Value.real() ;
  				Input[j + NewSamples * i][1] = Value.imag() ; 
  			}
  		}
  	}
}

void SpectralBIOs::AssembleBIO(map < pair< int , int > , complex < double > > * BIO_M) {
  	int Index1 ;
  	int Index2 ;
  	int IndexT ;
	for(int i = -nModes; i<=nModes; i++) {
		Index1 = GetIndex(i,Samples) ;
		for(int j = -nModes; j<=nModes; j++) {
			Index2 = GetIndex(j,Samples) ;
			IndexT = Index2 + Index1 * NewSamples ;
			(*BIO_M)[make_pair(i,-j)] = (*BIO_M)[make_pair(i,-j)] + 1.0/(double)(NewSamples * NewSamples)*complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
		}
	}
}

// === Single Layer Operator V ===
struct OperatorV {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	double WaveNumber ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	// === Laplace ===
	complex< double > ComputeG1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return complex < double > (-1.0 / (4.0 * M_PI)*log(pow(length(P1,P2)/(2.0*sin(M_PI *(x-y))),2.0)),0.0) ;
	}
	//
	complex < double > ComputeG2(double x) {
		ParamBoundRep->GetBoundaryRepresentationDerivative(P1,x) ;
		return complex < double > (-1.0 / (4.0 * M_PI)*log(pow(length(P1)/(2.0 * M_PI),2.0)),0.0) ;
	}
	// === Helmholtz ===
	complex <double> ComputeM1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return complex<double>(0.0,0.25) * boost::math::cyl_hankel_1(0,WaveNumber*length(P1,P2)) + 1.0/(2.0*M_PI)*log(length(P1,P2)) ;
	}
	complex <double> ComputeM2(double x) {
		return complex<double>(0.0,0.25) - 1.0/(2.0 * M_PI) * (log(WaveNumber/2.0) + (boost::math::constants::euler<double>())) ;
	}
	// === Helmholtz New ===
	complex < double > ComputeM1_1(double x, double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		return -1.0/(4.0*M_PI)*(boost::math::cyl_bessel_j(0,WaveNumber*length(P1,P2))-1.0);
	}
	complex < double > ComputeM1_2(double x, double y) {
		return ComputeM1(x,y) - log(pow((2.0*sin(M_PI *(x-y))),2.0))*ComputeM1_1(x,y);
	}
	complex < double > ComputeM2_1(double x) {
		return complex < double > (0.0,0.0);
	}
	complex < double > ComputeM2_2(double x) {
		return ComputeM2(x);
	}
} ;

complex < double > GetG1(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeG1(x,y) ;
}

complex < double > GetG2(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeG2(x) ;
}

complex < double > GetM1(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM1(x,y) ;
}

complex < double > GetM2(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM2(x) ;
}

void SpectralBIOs::BuildSpectralV() {
	V.clear();
	// 
	OperatorV OpV ;
	OpV.ParamBoundRep = PBR ;
  	//
	SampleBIO(GetG1,GetG2,&OpV);
 	fftw_execute(PlanSpectralBIOs) ;
	AssembleBIO(&V);
	//
  	for(int i = -nModes; i<=nModes; i++) {
  		if(i != 0) {
  			V[make_pair(i,i)] = V[make_pair(i,i)] + 1.0 / (4.0 * M_PI * abs((double)(i))) ;
  		}
  	}
	//
  	if(CaseWaveNumber=="Helmholtz") {
		OpV.WaveNumber = WaveNumber ;
		SampleBIO(GetM1,GetM2,&OpV);
 		fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&V);
  	}
	cout << V[make_pair(0,0)] << "\n";
}

// === Helmholtz New ===
complex < double > ComputeM1_1(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM1_1(x,y) ;
}
complex < double > ComputeM1_2(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM1_2(x,y) ;
}
complex < double > ComputeM2_1(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM2_1(x) ;
}
complex < double > ComputeM2_2(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	return OV->ComputeM2_2(x) ;
}

double CoeffLog(int x) {
	double v = 0;
	if (x!=0) {
		v = -1.0/abs((double)(x)); 
	}
	return v;
}

void SpectralBIOs::BuildSpectralV2() {
	//
	V.clear();
	// 
	OperatorV OpV ;
	OpV.ParamBoundRep = PBR ;
  	//
	SampleBIO(GetG1,GetG2,&OpV);
 	fftw_execute(PlanSpectralBIOs) ;
	AssembleBIO(&V);
	//
  	for(int i = -nModes; i<=nModes; i++) {
  		if(i != 0) {
  			V[make_pair(i,i)] = V[make_pair(i,i)] + 1.0 / (4.0 * M_PI * abs((double)(i))) ;
  		}
  	}
	//
	if(CaseWaveNumber=="Helmholtz") { 	
		OpV.WaveNumber = WaveNumber ;
		map < pair< int , int > , complex < double > > Aux;
		SampleBIO(ComputeM1_1,ComputeM2_1,&OpV);
		fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&Aux);
		//
		complex < double > Value = complex < double > (0.0,0.0);
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
				Value = complex < double > (0.0,0.0);
				for(int k = -Samples; k<=Samples; k++) {
					//if (i-j-k<=nModes && i-j-k>=nModes){
					Value = Value + CoeffLog(k+j) * Aux[make_pair(i-j-k,-k)] ;
					//}
				}
				V[make_pair(i,j)] = V[make_pair(i,j)] + Value;
  			}
  		}
		//
		SampleBIO(ComputeM1_2,ComputeM2_2,&OpV);
	 	fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&V);
	}
	cout << V[make_pair(0,0)] << "\n";
}

// === Double Layer Operator K ===
struct OperatorK {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	double WaveNumber ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	Point Normal ;
	double NtV = 0;
	//
	bool method;
	// === Laplace ===
	complex < double > ComputeK1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		ParamBoundRep->GetNormalVector(Normal,y) ;
		NtV = ParamBoundRep->NormOfTangentVector(y) ;
		//complex < double > Value = complex < double > (0.0,0.0);
		return complex < double > (1.0/(2.0 * M_PI) * (((P1-P2) & Normal))/(pow(length(P1,P2),2.0))*NtV,0.0);
	}
	//
	complex < double > ComputeK2(double x) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentationDerivative(P2, x) ;
		ParamBoundRep->GetBoundaryRepresentationSecondDerivative(P3, x) ;
		ParamBoundRep->GetNormalVector(Normal,x) ;
		NtV = ParamBoundRep->NormOfTangentVector(x ) ;
		return complex < double > (1.0/(4.0 * M_PI) * ((Normal & P3))/(pow(length(P2),2.0))*NtV,0.0);
	}
	// === Helmholtz ===
	complex< double > ComputeH1(double x, double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		ParamBoundRep->GetNormalVector(Normal,y) ;
		NtV = ParamBoundRep->NormOfTangentVector(y) ;
		return complex< double >(0.0,0.25*WaveNumber*NtV)*boost::math::cyl_hankel_1(1,WaveNumber*length(P1,P2))*((P1-P2)&Normal)/(length(P1,P2));
	}
	//
	complex < double > ComputeH2(double x) {
		return ComputeK2(x);
	}
	// === Helmholtz New ===
	complex < double > ComputeH1_1(double x, double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		return -0.5*WaveNumber*length(P1,P2)*boost::math::cyl_bessel_j(1,WaveNumber*length(P1,P2))*ComputeK1(x,y);
	}
	complex < double > ComputeH1_2(double x, double y) {
		return ComputeH1(x,y) - log(pow((2.0*sin(M_PI *(x-y))),2.0))*ComputeH1_1(x,y);
	}
	//
	complex < double > ComputeH2_1(double x) {
		return complex < double > (0.0,0.0);
	}
	complex < double > ComputeH2_2(double x) {
		return ComputeH2(x);
	}
} ;

// === Laplace ===
complex < double > GetK1(double x, double y, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeK1(x,y) ;
}
complex < double > GetK2(double x, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeK2(x) ;
}
// === Helmholtz ===
complex < double > GetH1(double x, double y, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH1(x,y) ;
}
complex < double > GetH2(double x, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH2(x) ;
}

void SpectralBIOs::BuildSpectralK() {
	K.clear();
	// 
	OperatorK OpK ;
	OpK.ParamBoundRep = PBR ;
	OpK.WaveNumber = WaveNumber ;
	//
	map < pair< int , int > , complex < double > > Aux  ;
	//
	if(CaseWaveNumber=="Laplace") { 
		SampleBIO(GetK1,GetK2,&OpK);
	 	fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&K);
	}
	//
	if(CaseWaveNumber=="Helmholtz") { 	
		SampleBIO(GetH1,GetH2,&OpK);
	 	fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&K);
	}
}

// === Helmholtz New ===
complex < double > ComputeH1_1(double x, double y, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH1_1(x,y) ;
}
complex < double > ComputeH1_2(double x, double y, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH1_2(x,y) ;
}
complex < double > ComputeH2_1(double x, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH2_1(x) ;
}
complex < double > ComputeH2_2(double x, void * Input) {
	OperatorK * OK = (OperatorK *) Input ;
	return OK->ComputeH2_2(x) ;
}

void SpectralBIOs::BuildSpectralK2() {
	K.clear();
	// 
	OperatorK OpK ;
	OpK.ParamBoundRep = PBR;
	//
	if(CaseWaveNumber=="Laplace") { 
		SampleBIO(GetK1,GetK2,&OpK);
		fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&K);
	}
	//
	if(CaseWaveNumber=="Helmholtz") { 	
		OpK.WaveNumber = WaveNumber ;
		//
		map < pair< int , int > , complex < double > > Aux;
		SampleBIO(ComputeH1_1,ComputeH2_1,&OpK);
		fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&Aux);
		//
		complex < double > Value = complex < double > (0.0,0.0);
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
				Value = complex < double > (0.0,0.0);
				for(int k = -Samples; k<=Samples; k++) {
					Value = Value + CoeffLog(k+j) * Aux[make_pair(i-j-k,-k)] ;
				}
				K[make_pair(i,j)] = Value;//K[make_pair(i,j)] + Value;
  			}
  		}
		//
		SampleBIO(ComputeH1_2,ComputeH2_2,&OpK);
	 	fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&K);
	}
	// cout << setprecision(10) << K[make_pair(5,5)] << "\n";
}

// === Adjoint Double Layer K' ===
void SpectralBIOs::BuildSpectralKt() {
  	// 
  	int Index1 ;
  	int Index2 ;
  	int IndexT ;
  	// Works for Laplace and Helmholtz
  	for(int i = -nModes; i<=nModes; i++) {
  		Index1 = GetIndex(i,Samples) ;
  		for(int j = -nModes; j<=nModes; j++) {
  			Kt[make_pair(i,j)] = K[make_pair(-j,-i)] ;	
  		}
  	}
}

struct TangentVector {
	ParametricBoundaryRepresentation * PBR ;
	Point P ;
	complex<double> Wave(double t) {	
		PBR->GetBoundaryRepresentationDerivative(P,t) ;
		return complex < double > (P[0],P[1]);
	}
} ;

complex < double > GetNu1(double t, void * Input) {
	TangentVector * TV = (TangentVector *) Input ;
	return complex < double > ((TV->Wave(t)).real(),0.0) ;
}

complex < double > GetNu2(double t, void * Input) {
	TangentVector * TV = (TangentVector *) Input ;
	return complex < double > ((TV->Wave(t)).imag(),0.0) ;
}

// === Hypersingular Operator W ===
void SpectralBIOs::BuildSpectralW() {
	W.clear();
	if(CaseWaveNumber=="Laplace") {
		// First Contribution
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  				}
  		}
		cout << setprecision(10) << W[make_pair(1,1)] << "\n";
	}
	//
	if(CaseWaveNumber=="Helmholtz") {
		// First Contribution
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  			}
  		}
		// Second Contribution
		TangentVector TV;
		TV.PBR = PBR;
		SpectralQuantity Nu1(nModes,Samples) ;
		SpectralQuantity Nu2(nModes,Samples) ;
		Nu1.ConvertToSpectral(GetNu1,&TV) ;
		Nu2.ConvertToSpectral(GetNu2,&TV) ;
		complex < double > Value = complex < double > (0.0,0.0);
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
				//
				Value = complex < double > (0.0,0.0);
				for(int k= -nModes; k<=nModes; k++) {
					for(int l=-nModes; l<=nModes; l++) {
						Value = Value + Nu1[k]*V[make_pair(-k+i,l+j)]*Nu1[l];
						Value = Value + Nu2[k]*V[make_pair(-k+i,l+j)]*Nu2[l];
					}		
				}
				W[make_pair(i,j)] = W[make_pair(i,j)] + pow(WaveNumber,2.0)*Value;
  			}
  		}
		cout << setprecision(10) << W[make_pair(-1,-1)] << "\n";
	}
	cout << setprecision(10) << W[make_pair(-1,-1)] << "\n";
}

// Log
complex < double > ComputeW1_1(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	Point N1;
	Point N2;
	OV->ParamBoundRep->GetNormalVector(N1,x) ;
	OV->ParamBoundRep->GetNormalVector(N2,y) ;
	return (N1&N2)*(-1.0/(4.0*M_PI) + (OV->ComputeM1_1(x,y))) ;
}
complex < double > ComputeW1_2(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	Point N;
	OV->ParamBoundRep->GetNormalVector(N,x) ;
	return   (N&N)*(-1.0/(4.0*M_PI) + (OV->ComputeM2_1(x))) ;
}
// Smooth
complex < double > ComputeW2_1(double x, double y, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	Point N1;
	Point N2;
	OV->ParamBoundRep->GetNormalVector(N1,x) ;
	OV->ParamBoundRep->GetNormalVector(N2,y) ;
	return (N1&N2)*((OV->ComputeM1_2(x,y)) + (OV->ComputeG1(x,y))) ;
}
complex < double > ComputeW2_2(double x, void * Input) {
	OperatorV * OV = (OperatorV *) Input ;
	Point N;
	OV->ParamBoundRep->GetNormalVector(N,x) ;
	return   (N&N)*((OV->ComputeM2_2(x))   + (OV->ComputeG2(x ))) ;
}

// === Hypersingular Operator W ===
void SpectralBIOs::BuildSpectralW2() {
	W.clear();
	if(CaseWaveNumber=="Laplace") {
		// First Contribution
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  			}
  		}
	}
	//
	if(CaseWaveNumber=="Helmholtz") {
		// First Contribution
		/*
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  			}
  		}
		*/
		// Second Contribution: 
		OperatorV OpV ;
		OpV.ParamBoundRep = PBR ;
		OpV.WaveNumber = WaveNumber ;
		//
		map < pair< int , int > , complex < double > > Aux1;
		map < pair< int , int > , complex < double > > Aux2;
		// Log
		SampleBIO(ComputeW1_1,ComputeW1_2,&OpV);
		fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&Aux1);
		// Smooth
		SampleBIO(ComputeW2_1,ComputeW2_2,&OpV);
	 	fftw_execute(PlanSpectralBIOs) ;
		AssembleBIO(&Aux2);
		//
		complex < double > Value = complex < double > (0.0,0.0);
		//
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
				Value = complex < double > (0.0,0.0);
				for(int k = -nModes; k<=nModes; k++) {
					if (i-j-k<=nModes && i-j-k>=nModes){
						Value = Value + CoeffLog(k+j) * Aux1[make_pair(i-j-k,-k)] ;
					}
				}
				W[make_pair(i,j)] = -pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] + pow(WaveNumber,2.0)*(Value+Aux2[make_pair(i,j)]);
  			}
  		}
	}
	cout << setprecision(10) << W[make_pair(-1,-1)] << "\n";
}
// === Mass Matrix === Just the Identity!
void SpectralBIOs::BuildMassMatrix() {
	for(int i = -nModes; i<=nModes; i++) {
		for(int j = -nModes; j<=nModes; j++) {
			M[make_pair(i,j)] = complex < double > (0.0,0.0);
			if (i==j) {
				M[make_pair(i,j)] = complex < double > (1.0,0.0);
			}
		}
	}
}

struct OperatorT {
	// T_{1,2}
	ParametricBoundaryRepresentation * ParamBoundRep1 ;
	ParametricBoundaryRepresentation * ParamBoundRep2 ;
	Point PX ;
	Point PY ;
	Point NX ;
	Point NY ;
	double JrX ;
	double JrY ;
	string CaseWaveNumber ;
	double WaveNumber ;
	//
	complex < double > Value ;
	complex < double > DirichletTraceDoubleLayer(double x, double y) {
		ParamBoundRep1->GetBoundaryRepresentation(PX, x) ;
		ParamBoundRep2->GetBoundaryRepresentation(PY, y) ;
		ParamBoundRep2->GetNormalVector(NY,y) ;
		JrY = (ParamBoundRep2->NormOfTangentVector(y)) ;
		if (CaseWaveNumber=="Laplace") { 
			return complex < double > (1.0/(2.0*M_PI) * (((PX-PY)&NY))*pow(length(PX,PY),-2.0) * JrY,0.0) ;
		} else if (CaseWaveNumber=="Helmholtz") { 
			return 0.0;
		} 
		return Value;
	}
	complex < double > DirichletTraceSingleLayer(double x, double y) {
		ParamBoundRep1->GetBoundaryRepresentation(PX, x) ;
		ParamBoundRep2->GetBoundaryRepresentation(PY, y) ;
		if(CaseWaveNumber=="Laplace") { 
			return -1.0/(2.0*M_PI) * log(length(PX,PY)); 
		} else if (CaseWaveNumber=="Helmholtz") { 
			return complex < double > (0.0,0.25)*boost::math::cyl_hankel_1(0,WaveNumber*length(PX,PY)) ;
		} else {
			assert(false) ;
		}
	}
	complex < double > NeumannTraceDoubleLayer(double x, double y) {
		ParamBoundRep1->GetBoundaryRepresentation(PX, x) ;
		ParamBoundRep2->GetBoundaryRepresentation(PY, y) ;
		ParamBoundRep1->GetNormalVector(NX,x) ;
		ParamBoundRep2->GetNormalVector(NY,y) ;
		JrX = (ParamBoundRep1->NormOfTangentVector(x)) ;
		JrY = (ParamBoundRep2->NormOfTangentVector(y)) ;
		if (CaseWaveNumber=="Laplace") { 
			return  -1.0/(2.0*M_PI) * (-(NX&NY)*pow(length(PX,PY),-2.0) + 2.0*((((PX-PY)&NX))*(((PX-PY)&NY)))*pow(length(PX,PY),-4.0)) * JrX * JrY ;
		} else if (CaseWaveNumber=="Helmholtz") { 
			return 0.0 ;
		} else {
			assert(false) ;
		}
	}
	complex < double > NeumannTraceSingleLayer(double x, double y) {
		ParamBoundRep1->GetBoundaryRepresentation(PX, x) ;
		ParamBoundRep2->GetBoundaryRepresentation(PY, y) ;
		ParamBoundRep1->GetNormalVector(NX,x) ;
		JrX = (ParamBoundRep2->NormOfTangentVector(x)) ;
		return  -1.0/(2.0*M_PI) * (((PX-PY)&NX))*pow(length(PX,PY),-2.0) * JrX ;
	}
} ;

void SpectralBIOs::SetBoundaryCrossInteraction(ParametricBoundaryRepresentation * PBR_) {
	PBR_C = PBR_ ;
}

void SpectralBIOs::BuildCrossInteraction() {
	//
	OperatorT OpT ;
	OpT.ParamBoundRep1 = PBR   ;
	OpT.ParamBoundRep2 = PBR_C ;
  	//
    double I1 ; 
  	double I2 ;
    // 
  	int Index1 ;
  	int Index2 ;
  	int IndexT ;
  	//
  	int Pos1 ;
  	int Pos2 ;
  	// Values
  	complex < double> Value = 0.0 ;
  	complex < double > ComplexValue = complex <double>(0.0,0.0) ;
  	// Par of Int
  	pair <int , int > PairOfInt ;
  	// Aux Map
  	map < pair< int , int > , complex < double > > Aux ;
  	for(int B = 0; B<4; B++) {
  		for(int i = 0; i<NewSamples; i++) {
  			I1 = (double)(i)/(double)(NewSamples) ;
			for(int j = 0; j<NewSamples; j++) {
				I2 = (double)(j)/(double)(NewSamples) ;
				if(B==0) {
					Value = OpT.DirichletTraceDoubleLayer(I1,I2) ;
				}
				else if (B==1) {
					Value = OpT.DirichletTraceSingleLayer(I1,I2) ;
				}
				else if (B==2) {
					Value = OpT.NeumannTraceDoubleLayer(I1,I2) ;
				}
				else {
					Value = OpT.NeumannTraceSingleLayer(I1,I2) ;
				}
				Input[j + NewSamples * i][0] = Value.real() ;
  				Input[j + NewSamples * i][1] = Value.imag(); 
  			}
  		}
  		//
 		fftw_execute(PlanSpectralBIOs) ;
		//
  		for(int i = -nModes; i<=nModes; i++) {
  			Index1 = GetIndex(i,Samples) ;
  			for(int j = -nModes; j<=nModes; j++) {
  					Index2 = GetIndex(j,Samples) ;
  					IndexT = Index2 + Index1 * NewSamples ;
  					Aux[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  			}
  		}
  		T[B] = Aux ;
	}
}

map < pair< int , int > , complex < double > > * SpectralBIOs::GetV() {
	return &V ;
}

map < pair< int , int > , complex < double > > * SpectralBIOs::GetW() {
	return &W ; 
}

map < pair< int , int > , complex < double > > * SpectralBIOs::GetK() {
	return &K;
}

map < pair< int , int > , complex < double > > * SpectralBIOs::GetKt() {
	return &Kt ; 
}

map < pair< int , int > , complex < double > > * SpectralBIOs::GetT(int B) {
	return &(T[B]) ; 
}

SpectralBIOs::~SpectralBIOs() {
	fftw_destroy_plan(PlanSpectralBIOs) ;
	fftw_free(Input)  ;
  	fftw_free(Output) ;
  	fftw_destroy_plan(PlanSpectralBIOs_1D) ;
	fftw_free(Input_1D)  ;
  	fftw_free(Output_1D) ;
}


