#ifndef NEUMANN_PROBLEM_HPP
#define NEUMANN_PROBLEM_HPP

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

class NeumannProblem {
private:
	SpectralBIOs * BIOs ;
	MatrixXcd W ;
	VectorXcd f ;
	VectorXcd x ;
	
public:
	NeumannProblem(SpectralBIOs & BIOs_) ;
	void BuildMatrix() ;
	void BuildRHSDirectMethod(SpectralQuantity & SQ) ;
	void BuildRHSIndirectMethod(SpectralQuantity & SQ) ;
	void Solve() ;
	void GetSolution(SpectralQuantity & SQ) ;
	
} ;
#endif

				