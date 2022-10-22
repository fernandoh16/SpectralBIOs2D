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
	BuildSpectralKt() ;
	BuildSpectralW() ;
	//BuildCrossInteraction() ;
}

void SpectralBIOs::Compute_FFT_SpectralBIOS(map < pair< int , int > , complex < double > > * M){
	// 
	OperatorV OpV ;
	OpV.ParamBoundRep = PBR ;
  	//
    double I1 ; 
  	double I2 ;
  	//
  	for(int i = 0; i<NewSamples; i++) {
  		I1 = (double)(i)/(double)(NewSamples) ;
		for(int j = 0; j<NewSamples; j++) {
			I2 = (double)(j)/(double)(NewSamples) ;
			if(i != j) {
				Input[j + NewSamples * i][0] = OpV.ComputeG1(I1,I2) ;
  				Input[j + NewSamples * i][1] = 0.0 ; 
  			}
  			if (i == j) {
  				Input[j + NewSamples * i][0] = OpV.ComputeG2(I1) ;
  				Input[j + NewSamples * i][1] = 0.0 ; 
  			}
  		}
  	}
 	fftw_execute(PlanSpectralBIOs) ;
  	// 
  	int Index1 ;
  	int Index2 ;
  	int IndexT ;
  	//
  	for(int i = -nModes; i<=nModes; i++) {
  		Index1 = GetIndex(i,Samples) ;
  		for(int j = -nModes; j<=nModes; j++) {
  				Index2 = GetIndex(j,Samples) ;
  				IndexT = Index2 + Index1 * NewSamples ;
  				V[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  		}
  	}
 	//
	
  	for(int i = -nModes; i<=nModes; i++) {
  		if(i != 0) {
  			V[make_pair(i,i)] = V[make_pair(i,i)] + 1.0 / (4.0 * M_PI * abs((double)(i))) ;
  		}
  	}	
	
}

// === Single Layer Operator V ===
// Laplace - Checked
/*
struct OperatorV {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	// 
	double ComputeG1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return -1.0 / (4.0 * M_PI) * log(pow(length(P1,P2)/(2.0*sin(M_PI *(x-y))),2.0)) ;
	}
	//
	double ComputeG2(double x) {
		ParamBoundRep->GetBoundaryRepresentationDerivative(P1,x) ;
		return -1.0 / (4.0 * M_PI) * log(pow(length(P1)/(2.0 * M_PI),2.0)) ;
	}
} ;
// Helmholtz
struct OperatorV_Helm {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	//
	double Wavenumber ;
	// 
	// === Helmholtz == 
	complex <double> ComputeG1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return complex<double>(0.0,0.25) * boost::math::cyl_hankel_1(0,Wavenumber*length(P1,P2)) + 1.0/(2.0*M_PI)*log(length(P1,P2)) ;
	}
	complex <double> ComputeG2(double x) {
		return complex<double>(0.0,0.25) - 1.0/(2.0 * M_PI) * (log(Wavenumber/2.0) + (boost::math::constants::euler<double>())) ;
	}
} ;
*/

// === Single Layer Operator V ===
struct OperatorV {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	double WaveNumber ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	// === Laplace ===
	double ComputeG1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return -1.0 / (4.0 * M_PI) * log(pow(length(P1,P2)/(2.0*sin(M_PI *(x-y))),2.0)) ;
	}
	//
	double ComputeG2(double x) {
		ParamBoundRep->GetBoundaryRepresentationDerivative(P1,x) ;
		return -1.0 / (4.0 * M_PI) * log(pow(length(P1)/(2.0 * M_PI),2.0)) ;
	}
	// === Helmholtz == 
	complex <double> ComputeM1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1,x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2,y) ;
		return complex<double>(0.0,0.25) * boost::math::cyl_hankel_1(0,WaveNumber*length(P1,P2)) + 1.0/(2.0*M_PI)*log(length(P1,P2)) ;
	}
	complex <double> ComputeM2(double x) {
		return complex<double>(0.0,0.25) - 1.0/(2.0 * M_PI) * (log(WaveNumber/2.0) + (boost::math::constants::euler<double>())) ;
	}
} ;

void SpectralBIOs::BuildSpectralV() {
	// 
	OperatorV OpV ;
	OpV.ParamBoundRep = PBR ;
  	//
    double I1 ; 
  	double I2 ;
  	//
  	for(int i = 0; i<NewSamples; i++) {
  		I1 = (double)(i)/(double)(NewSamples) ;
		for(int j = 0; j<NewSamples; j++) {
			I2 = (double)(j)/(double)(NewSamples) ;
			if(i != j) {
				Input[j + NewSamples * i][0] = OpV.ComputeG1(I1,I2) ;
  				Input[j + NewSamples * i][1] = 0.0 ; 
  			}
  			if (i == j) {
  				Input[j + NewSamples * i][0] = OpV.ComputeG2(I1) ;
  				Input[j + NewSamples * i][1] = 0.0 ; 
  			}
  		}
  	}
 	fftw_execute(PlanSpectralBIOs) ;
  	// 
  	int Index1 ;
  	int Index2 ;
  	int IndexT ;
  	//
  	for(int i = -nModes; i<=nModes; i++) {
  		Index1 = GetIndex(i,Samples) ;
  		for(int j = -nModes; j<=nModes; j++) {
  				Index2 = GetIndex(j,Samples) ;
  				IndexT = Index2 + Index1 * NewSamples ;
  				V[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  		}
  	}
 	//
	
  	for(int i = -nModes; i<=nModes; i++) {
  		if(i != 0) {
  			V[make_pair(i,i)] = V[make_pair(i,i)] + 1.0 / (4.0 * M_PI * abs((double)(i))) ;
  		}
  	}
	
	// Sanity Test 
	/*
	int N_Quad = 1000;
	int m = -3;
	int n = 2;
	complex < double > Value = complex < double > (0.0,0.0);
	double N_Q1, N_Q2;
	double dh = 1.0/ (double)(N_Quad);
	for(int i = 0; i<N_Quad; i++) {	
		N_Q1 = (double)(i)/(double)(N_Quad);
		for(int j = 0; j<N_Quad; j++) {
			N_Q2 = (double)(j)/(double)(N_Quad);
			if(i!=j) {
				Value = Value + pow(dh,2)*OpV.ComputeG1(N_Q1,N_Q2) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
			}
			if(i==j){
				Value = Value + pow(dh,2)*OpV.ComputeG2(N_Q1) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
			}
		}
	}
	cout << Value << " " << V[make_pair(m,n)];
	*/
	//
  	if(CaseWaveNumber=="Helmholtz") {
		// OperatorV_Helm OpV_Helm ;
		// OpV_Helm.ParamBoundRep = PBR ;
		OpV.WaveNumber = WaveNumber ;
  		//
    	// double I1 ; 
  		// double I2 ;
  		//
  		for(int i = 0; i<NewSamples; i++) {
  			I1 = (double)(i)/(double)(NewSamples) ;
			for(int j = 0; j<NewSamples; j++) {
				I2 = (double)(j)/(double)(NewSamples) ;
				if(i != j) {
					Input[j + NewSamples * i][0] = (OpV.ComputeM1(I1,I2)).real() ;
  					Input[j + NewSamples * i][1] = (OpV.ComputeM1(I1,I2)).imag() ; 
  				}
  				if (i == j) {
  					Input[j + NewSamples * i][0] = (OpV.ComputeM2(I1)).real() ;
  					Input[j + NewSamples * i][1] = (OpV.ComputeM2(I1)).imag() ; 
  				}
  			}
  		}
 		fftw_execute(PlanSpectralBIOs) ;
  		// 
  		// int Index1 ;
  		// int Index2 ;
  		// int IndexT ;
  		//
  		for(int i = -nModes; i<=nModes; i++) {
  			Index1 = GetIndex(i,Samples) ;
  			for(int j = -nModes; j<=nModes; j++) {
  				Index2 = GetIndex(j,Samples) ;
  				IndexT = Index2 + Index1 * NewSamples ;
  				V[make_pair(i,-j)] = V[make_pair(i,-j)] + 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
				//V[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  			}
  		}
		// Sanity Test (Only for the Helmholtz contribution)
		/*
		int N_Quad = 5000;
		int m = 50;
		int n = 50;
		complex < double > Value = complex < double > (0.0,0.0);
		double N_Q1, N_Q2;
		double dh = 1.0/ (double)(N_Quad);
		cout << OpV_Helm.Wavenumber << "\n";
		for(int i = 0; i<N_Quad; i++) {	
			N_Q1 = (double)(i)/(double)(N_Quad);
			for(int j = 0; j<N_Quad; j++) {
				N_Q2 = (double)(j)/(double)(N_Quad);
				if(i!=j) {
					Value = Value + pow(dh,2)*OpV_Helm.ComputeG1(N_Q1,N_Q2) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
				}
				if(i==j){
					Value = Value + pow(dh,2)*OpV_Helm.ComputeG2(N_Q1) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
				}
			}
		}
		cout << Value << " " << V[make_pair(m,n)];
		*/
  	}
}

// === Double Layer Operator K ===
/*
struct OperatorK {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	Point Normal ;
	// 
	double ComputeK1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		ParamBoundRep->GetNormalVector(Normal,y) ;
		return 1.0/(2.0 * M_PI) * (((P1-P2) & Normal)) /(pow(length(P1,P2),2.0)) * (ParamBoundRep->NormOfTangentVector(y)) ;
	}
	//
	double ComputeK2(double x) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentationDerivative(P2, x) ;
		ParamBoundRep->GetBoundaryRepresentationSecondDerivative(P3, x) ;
		ParamBoundRep->GetNormalVector(Normal,x) ;
		return 1.0/(4.0 * M_PI) * ((Normal & P3)) /(pow(length(P2),2.0)) * (ParamBoundRep->NormOfTangentVector(x));
	}
} ;
*/
// === Double Layer Operator K ===

struct OperatorK {
	ParametricBoundaryRepresentation * ParamBoundRep ;
	double WaveNumber ;
	Point P1 ;
	Point P2 ;
	Point P3 ;
	Point Normal ;
	// === Laplace ===
	double ComputeK1(double x , double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		ParamBoundRep->GetNormalVector(Normal,y) ;
		return 1.0/(2.0 * M_PI) * (((P1-P2) & Normal))/(pow(length(P1,P2),2.0)) ;
	}
	//
	double ComputeK2(double x) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentationDerivative(P2, x) ;
		ParamBoundRep->GetBoundaryRepresentationSecondDerivative(P3, x) ;
		ParamBoundRep->GetNormalVector(Normal,x) ;
		return 1.0/(4.0 * M_PI) * ((Normal & P3))/(pow(length(P2),2.0));
	}
	// === Helmholtz == 
	complex<double> ComputeH1(double x, double y) {
		ParamBoundRep->GetBoundaryRepresentation(P1, x) ;
		ParamBoundRep->GetBoundaryRepresentation(P2, y) ;
		ParamBoundRep->GetNormalVector(Normal,y) ;
		complex <double> Value ;
		Value = complex<double>(0.0,0.25*WaveNumber)*boost::math::cyl_hankel_1(1,WaveNumber*length(P1,P2))*((P1-P2)&Normal)/(length(P1,P2));
		return  Value - ComputeK1(x,y) ;
	}
	//
	complex<double> ComputeH2(double x) {
		return complex < double > (0.0,0.0) ;
	}
} ;

void SpectralBIOs::BuildSpectralK() {
	// 
	OperatorK OpK ;
	OpK.ParamBoundRep = PBR ;
	//
    double I1 ; 
  	double I2 ;
  	//
	int Index1 ;
	int Index2 ;
	int IndexT ;
	//
	if(CaseWaveNumber=="Laplace") { 
		for(int i = 0; i<NewSamples; i++) {
  			I1 = (double)(i)/(double)(NewSamples) ;
			for(int j = 0; j<NewSamples; j++) {
				I2 = (double)(j)/(double)(NewSamples) ;
				if(i != j) {
					Input[j + NewSamples * i][0] = OpK.ComputeK1(I1,I2) ;
  					Input[j + NewSamples * i][1] = 0.0 ;
  				}
  				if (i == j) {
  					Input[j + NewSamples * i][0] = OpK.ComputeK2(I1) ;
  					Input[j + NewSamples * i][1] = 0.0 ;
  				}
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
  				K[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  			}
  		}
	}
	//
	if(CaseWaveNumber=="Helmholtz") { 	
		//
		OpK.WaveNumber = WaveNumber ;
		for(int i = 0; i<NewSamples; i++) {
  			I1 = (double)(i)/(double)(NewSamples) ;
			for(int j = 0; j<NewSamples; j++) {
				I2 = (double)(j)/(double)(NewSamples) ;
				if(i != j) {
					Input[j + NewSamples * i][0] = (OpK.ComputeH1(I1,I2)).real() ;
  					Input[j + NewSamples * i][1] = (OpK.ComputeH1(I1,I2)).imag() ;
  				}
  				if (i == j) {
  					Input[j + NewSamples * i][0] = (OpK.ComputeH2(I1)).real() ;
  					Input[j + NewSamples * i][1] = (OpK.ComputeH2(I1)).imag() ;
  				}
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
  				K[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  			}
  		}
		// Sanity Test (Only for the Helmholtz contribution)
		int N_Quad = 500;
		int m =  0;
		int n =  0;
		complex < double > Value = complex < double > (0.0,0.0);
		double N_Q1, N_Q2;
		double dh = 1.0/ (double)(N_Quad);
		// cout << OpV_Helm.Wavenumber << "\n";
		for(int i = 0; i<N_Quad; i++) {	
			N_Q1 = (double)(i)/(double)(N_Quad);
			for(int j = 0; j<N_Quad; j++) {
				N_Q2 = (double)(j)/(double)(N_Quad);
				if(i!=j) {
					Value = Value + pow(dh,2)*OpK.ComputeH1(N_Q1,N_Q2) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
				}
				if(i==j){
					Value = Value + pow(dh,2)*OpK.ComputeH2(N_Q1) * exp(complex < double > (0,2*M_PI*(double)(n)*N_Q1-2*M_PI*(double)(m)*N_Q2));
				}
			}
		}
		cout << Value << " " << K[make_pair(m,n)] << "\n";
	}
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

// === Hypersingular Operator W ===
void SpectralBIOs::BuildSpectralW() {
	// 
	if(CaseWaveNumber=="Laplace") {
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  			}
  		}
	}
	if(CaseWaveNumber=="Helmholtz") {
		// First Contribution
  		for(int i = -nModes; i<=nModes; i++) {
  			for(int j = -nModes; j<=nModes; j++) {
  				W[make_pair(i,j)] = pow(2.0 * M_PI, 2.0) * (double)(i) * (double)(j) * V[make_pair(i,j)] ;
  			}
  		}
		/*
		// Second Contribution
	    double I1 ; 
	  	double I2 ;
	  	//
		int Index1 ;
		int Index2 ;
		int IndexT ;
		//
		Point Normal1 ;
		Point Normal2 ;
		for(int i = 0; i<NewSamples; i++) {
  			I1 = (double)(i)/(double)(NewSamples) ;
			ParamBoundRep->GetNormalVector(Normal1,I1) ;
			for(int j = 0; j<NewSamples; j++) {
				I2 = (double)(j)/(double)(NewSamples) ;
				ParamBoundRep->GetNormalVector(Normal2,I2) ;
				Input[j + NewSamples * i][0] = (Normal1 & Normal2);
  				Input[j + NewSamples * i][1] = 0.0
  			}
  		}
		//
 		fftw_execute(PlanSpectralBIOs) ;
  		//
		map < pair< int , int > , complex < double > > M_Normal ;
  		for(int i = -nModes; i<=nModes; i++) {
  			Index1 = GetIndex(i,Samples) ;
  			for(int j = -nModes; j<=nModes; j++) {
  				Index2 = GetIndex(j,Samples) ;
  				IndexT = Index2 + Index1 * NewSamples ;
  				K[make_pair(i,-j)] = 1.0 / (double)(NewSamples * NewSamples) * complex<double>(Output[IndexT][0],Output[IndexT][1]) ;
  			}
  		}
		*/
	}
}

// === Mass Matrix ===
void SpectralBIOs::BuildMassMatrix() {
  	double S ;
  	double V ;
  	for (int i = 0; i<NewSamples; i++) {
  		S = (double)(i)/(double)(NewSamples) ;
  		V = 1./(PBR->NormOfTangentVector(S)) ;
    	Input_1D[i][0] =   V ;
    	Input_1D[i][1] = 0.0 ;
  	}
	// 
  	fftw_execute(PlanSpectralBIOs_1D) ;
  	// Store data in SQ
	map < int, complex < double > > SQ ;
  	int Index ;
  	for (int i = 0; i<=Samples; i++) {
  		if(i == 0) {
  			SQ[0]  = 1.0/(double)(NewSamples) * complex < double > (Output_1D[0][0],Output_1D[0][1]) ;
  		}
  		else {
			Index  = 2 * Samples - i + 1;
			SQ[-i] = 1.0/(double)(NewSamples) * complex < double > (Output_1D[Index][0],Output_1D[Index][1]) ;
  			Index  = i ;
			SQ[+i] = 1.0/(double)(NewSamples) * complex < double > (Output_1D[Index][0],Output_1D[Index][1]) ;
		}
  	}
	for(int i = -nModes; i<=nModes; i++) {
		for(int j = -nModes; j<=nModes; j++) {
			M[make_pair(i,j)] = SQ[i-j] ;
		}
	}
	//
	/*
	int N_Quad = 5000;
	int m =  1 ;
	int n = -5 ;
	complex < double > Value = complex < double > (0.0,0.0);
	double N_Q;
	double dh = 1.0/ (double)(N_Quad);
	for(int i = 0; i<N_Quad; i++) {	
		N_Q = (double)(i)/(double)(N_Quad);
		V = 1./(PBR->NormOfTangentVector(N_Q)) ;
		Value = Value + dh*V* exp(complex < double > (0,2*M_PI*(double)(n)*N_Q-2*M_PI*(double)(m)*N_Q));
	}
	cout << Value << " " << M[make_pair(m,n)] << "\n";
	*/
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
	complex < double > DirichletTraceDoubleLayer(double x, double y) {
		ParamBoundRep1->GetBoundaryRepresentation(PX, x) ;
		ParamBoundRep2->GetBoundaryRepresentation(PY, y) ;
		ParamBoundRep2->GetNormalVector(NY,y) ;
		JrY = (ParamBoundRep2->NormOfTangentVector(y)) ;
		if (CaseWaveNumber=="Laplace") { 
			return complex < double > (1.0/(2.0*M_PI) * (((PX-PY)&NY))*pow(length(PX,PY),-2.0) * JrY,0.0) ;
		} else if (CaseWaveNumber=="Helmholtz") { 
			return 0.0 ;
		} else {
			assert(false) ;
		}
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
		return  -1.0/(2.0*M_PI) * (-(NX&NY)*pow(length(PX,PY),-2.0) + 2.0*((((PX-PY)&NX))*(((PX-PY)&NY)))*pow(length(PX,PY),-4.0)) * JrX * JrY ;
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


