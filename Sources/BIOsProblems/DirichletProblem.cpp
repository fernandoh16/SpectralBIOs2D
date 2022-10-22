#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>
#include <map> 

// Sources
#include "DirichletProblem.hpp"
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"

// Eigen
#include <Eigen/Dense>
using namespace Eigen;

DirichletProblem::DirichletProblem(SpectralBIOs & BIOs_): 
	V(2*(BIOs_.GetNumberofModes())+1,2*(BIOs_.GetNumberofModes())+1), 
	f(2*(BIOs_.GetNumberofModes())+1), 
	x(2*(BIOs_.GetNumberofModes())+1) {
	BIOs = &BIOs_ ;
}

void DirichletProblem::BuildMatrix(){
	map < pair< int , int > , complex < double > > * Vmap ;
	Vmap = BIOs->GetV() ;
	int nModes = BIOs->GetNumberofModes() ;
	for(int i =-nModes; i<=nModes; i++) {
		for(int j =-nModes; j<=nModes; j++) {
			V(i + nModes,j + nModes) = (*Vmap)[make_pair(i,j)] ;
		}
	}	
}

void DirichletProblem::BuildRHSDirectMethod(SpectralQuantity & SQ1, SpectralQuantity & SQ2) {
	map < pair< int , int > , complex < double > > * Kmap ;
	Kmap = BIOs->GetK() ;
	int nModes = BIOs->GetNumberofModes() ;
	for(int i =-nModes; i<=nModes; i++) {
		f(i+nModes) = 0.5 * SQ1[i] ;
		for(int j = -nModes; j<=nModes; j++) {
			f(i+nModes) = f(i+nModes) + ((*Kmap)[make_pair(i,j)] * SQ2[j]) ;
		}
	}
	/*
	for(int i =-nModes; i<=nModes; i++) {
		for(int j =-nModes; j<=nModes; j++) {
			cout << (*Kmap)[make_pair(i,j)]  << "  ";
		}
		cout << "\n" ;
	}
	*/
}

void DirichletProblem::BuildRHSIndirectMethod(SpectralQuantity & SQ) {
	int nModes = BIOs->GetNumberofModes() ;
	for(int i = -nModes; i<=nModes; i++) {
		f(i+nModes) = SQ[i] ;
	}
}

void DirichletProblem::Solve() {
	// x = V.partialPivLU.solve(f);
	x = V.fullPivLu().solve(f);
}

void DirichletProblem::GetSolution(SpectralQuantity & SQ) {	
	int nModes = BIOs->GetNumberofModes() ;
	for(int i = -nModes; i<=nModes; i++) {
		SQ[i] = x(i+nModes) ;
	}
}


				