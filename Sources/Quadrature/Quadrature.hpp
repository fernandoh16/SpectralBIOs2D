#ifndef QUADRATUREHEADERDEF 
#define QUADRATUREHEADERDEF

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

using namespace std ;

class Quadrature {
private:
	string Dimension ;
	int NumberofPoints ;
	int NumberofCycles ;
	vector < double > QuadraturePoints  ;
	vector < double > QuadratureWeights ;
	
public:
	Quadrature() ;	
	void SetQuadrature(string Dimension_, int  NumberofPoints_, int NumberofCycles_) ;	
	double ComputeQuadrature(double (*p_func)(double x, void * Input), void * Input_, double x1 = -1.0, double x2 = +1.0) ;
	double ComputeQuadrature(double (*p_func)(double x, double y, void * Input), void * Input_, double x1 = -1.0, double x2 = +1.0, double y1 = -1.0, double y2 = +1.0) ;
	complex < double > ComputeQuadrature(complex < double > (*p_func)(double x, double y, void * Input), void * Input_, double x1, double x2, double y1, double y2) ;
	int GetNumberOfPoints() ;
	double GetQuadraturePoint(int i, double x1, double x2) ;
	double GetQuadratureWeight(int i, double x1, double x2) ;
} ; 
#endif

				