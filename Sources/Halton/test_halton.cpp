#include <cmath>
#include <string>
#include <iostream>
#include <cassert>
#include <ctime>

using namespace std ;

#include "halton.hpp"
#include "halton.cpp"

int main(int argc, char* argv[]) {
	int n_points = 10;
	int dim = 10;
	double * I;
	I = halton(10,dim);
	for (int i=0; i<dim; i++) {
		cout << I[i] << "\n";
	}
	cout << "Hola" << "\n";
	return 0 ;
}
