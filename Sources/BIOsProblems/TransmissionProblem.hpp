#ifndef TransmissionProblem_HPP
#define TransmissionProblem_HPP

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

class TransmissionProblem {
private:
	SpectralBIOs * BIOs ;
	MatrixXcd MTF ;
	VectorXcd f   ;
	VectorXcd x   ;
	
public:
	TransmissionProblem(SpectralBIOs & BIOs_) ;
	void BuildMatrix() ;
	void BuildRHS(SpectralQuantity & SQ1, SpectralQuantity & SQ1) ;
	void Solve() ;
	void GetSolution(SpectralQuantity & SQ) ;
	
} ;
#endif

				