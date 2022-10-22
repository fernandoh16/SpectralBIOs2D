#ifndef QoISOLUTION_HPP
#define QoISOLUTION_HPP

#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <iostream>

// Sources
#include "Point/Point.hpp"
#include "SpectralQuantity/SpectralQuantity.hpp"

using namespace std ;

class QoISolution {
private:
	int           NumberOfControlPoints ;
	vector <double>            ArcParam ;
	vector <complex <double> > Solution ;
	
public:
	QoISolution() ;
	QoISolution(int NumberOfControlPoints_) ;
	void ComputeSolution(SpectralQuantity & SQ) ;
	QoISolution & operator=(const QoISolution & QoI) ;
	QoISolution operator+(const QoISolution & QoI) const ;
	QoISolution operator-(const QoISolution & QoI) const ;
	QoISolution operator*(double a) const ;
	void PrintData() ;
	friend QoISolution operator*(double a, QoISolution const & QoI) ;
} ;
	QoISolution operator*(double a, QoISolution const & QoI) ;
#endif

				