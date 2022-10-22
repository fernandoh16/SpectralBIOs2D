#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>
#include <vector>

// Direct product
#include <direct_product.hpp>

// Sources
#include "BayesianInversion/BayesianInversion.hpp"

using namespace std ;

struct Measurements {
	double Value ; 
	Measurements(double Value_ = 0.0) {
		Value = Value_ ;
	}
	Measurements operator-(const Measurements & M) const {
		Measurements NewM ;
		NewM.Value = Value - M.Value ;
		return NewM ;
	}	
} ;

struct ForwardModel {
	double FM  ;
	Measurements M ;
	Measurements & operator()(vector < double > & y) {
		return M ;
	}
} ;

struct Potential {
	double operator()(const Measurements & M1, const Measurements & M2) const {
		Measurements NewM ;
		NewM = M1 - M2 ;
		return NewM.Value ;
	} 
} ;

int main(int argc, char* argv[]) {
	ForwardModel FM ;
	FM.M.Value = 0.0 ;
	Potential P ;
	Measurements M ;
	double A = 0.0 ;
	BayesianInversion <ForwardModel,Potential,Measurements,vector <double> >  BI(FM,P,M) ;
	vector < double > y ;
	y.push_back(0.0) ;
	auto out = BI(y) ;
	cout << get<2>(out) << "\n" ;
}
