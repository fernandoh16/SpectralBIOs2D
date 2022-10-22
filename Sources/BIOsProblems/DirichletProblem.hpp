#ifndef DIRICHLET_PROBLEM_HPP
#define DIRICHLET_PROBLEM_HPP

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

// Sources
#include "SpectralBIOs/SpectralBIOs.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"

// Eigen
#include <Eigen/Dense>
using namespace Eigen;

using namespace std ;

class DirichletProblem {
private:
	SpectralBIOs * BIOs ;
	MatrixXcd V ;
	VectorXcd f ;
	VectorXcd x ;
	
public:
	DirichletProblem(SpectralBIOs & BIOs_) ;
	void BuildMatrix() ;
	void BuildRHSDirectMethod(SpectralQuantity & SQ1, SpectralQuantity & SQ2) ;
	void BuildRHSIndirectMethod(SpectralQuantity & SQ) ;
	void Solve() ;
	void GetSolution(SpectralQuantity & SQ) ;
	
} ;
#endif

				