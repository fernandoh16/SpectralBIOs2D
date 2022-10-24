#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>
#include <map> 

// Sources
#include "NeumannProblem.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"

// Eigen
#include <Eigen/Dense>
using namespace Eigen;

NeumannProblem::NeumannProblem(SpectralBIOs & BIOs_): 
	W(2*(BIOs_.GetNumberofModes())+1,2*(BIOs_.GetNumberofModes())+1), 
	f(2*(BIOs_.GetNumberofModes())+1), 
	x(2*(BIOs_.GetNumberofModes())+1) {
	BIOs = &BIOs_ ;
}

void NeumannProblem::BuildMatrix(){
	map < pair< int , int > , complex < double > > * Wmap ;
	Wmap = BIOs->GetW() ;
	int nModes = BIOs->GetNumberofModes() ;
	for(int i =-nModes; i<=nModes; i++) {
		for(int j =-nModes; j<=nModes; j++) {
			W(i + nModes,j + nModes) = (*Wmap)[make_pair(i,j)] ;
		}
	}	
}

void NeumannProblem::BuildRHSDirectMethod(SpectralQuantity & SQ) {
	map < pair< int , int > , complex < double > > * Ktmap ;
	Ktmap = BIOs->GetKt() ;
	int nModes = BIOs->GetNumberofModes() ;
	for(int i =-nModes; i<=nModes; i++) {
		f(i+nModes) = 0.5 * SQ[i] ;
		for(int j = -nModes; j<=nModes; j++) {
			f(i+nModes) = f(i+nModes) - ((*Ktmap)[make_pair(i,j)] * SQ[j]) ;
		}
	}
}

void NeumannProblem::BuildRHSIndirectMethod(SpectralQuantity & SQ) {
	int nModes = BIOs->GetNumberofModes() ;
	for(int i = -nModes; i<=nModes; i++) {
		f(i+nModes) = SQ[i] ;
	}
}

void NeumannProblem::Solve() {
	// x = V.partialPivLU.solve(f);
	x = W.fullPivLu().solve(f);
	cout << (W*x-f).norm() << "\n";
}

void NeumannProblem::GetSolution(SpectralQuantity & SQ) {	
	int nModes = BIOs->GetNumberofModes() ;
	for(int i = -nModes; i<=nModes; i++) {
		SQ[i] = x(i+nModes) ;
	}
}


				