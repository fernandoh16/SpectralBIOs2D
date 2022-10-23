#ifndef KIDNEY_HEADERDEF 
#define KIDNEY_HEADERDEF
/*
#include <stdlib.h> 
#include <cmath> 

using namespace std ;
*/
void KidneyRep(Point & P, double t) {
	//
	double R = 0.2 ;
	//
	if (t>=0 && t<pow(6.0,-1.0)) {
		P[0] = 2*R+R*cos(6.0*M_PI*t) ;
		P[1] =     R*sin(6.0*M_PI*t) ;
	}
	if (t<pow(3.0,-1.0) && t>=pow(6.0,-1.0)) {
		P[0] = R*cos(6.0*M_PI*(t-pow(6.0,-1.0))) ;
		P[1] =-R*sin(6.0*M_PI*(t-pow(6.0,-1.0))) ;
	}
	if (t<pow(2.0,-1.0) && t>=pow(3.0,-1.0)) {
		P[0] =-2*R+R*cos(6.0*M_PI*(t-pow(3.0,-1.0))) ;
		P[1] =     R*sin(6.0*M_PI*(t-pow(3.0,-1.0))) ;
	}
	if (t<=1.0 && t>=pow(2.0,-1.0)) {
		P[0] =-3.0*R*cos(2.0*M_PI*(t-pow(2.0,-1.0))) ;
		P[1] =-3.0*R*sin(2.0*M_PI*(t-pow(2.0,-1.0))) ;
	}
}

void KidneyRepDer(Point & P, double t) {
	//
	double R = 0.2;
	//
	if (t>=0 && t<pow(6.0,-1.0)) {
		P[0] =-6.0*M_PI*R*sin(6.0*M_PI*t) ;
		P[1] = 6.0*M_PI*R*cos(6.0*M_PI*t) ;
	}
	if (t<pow(3.0,-1.0) && t>=pow(6.0,-1.0)) {
		P[0] =-6.0*M_PI*R*sin(6.0*M_PI*(t-pow(6.0,-1.0))) ;
		P[1] =-6.0*M_PI*R*cos(6.0*M_PI*(t-pow(6.0,-1.0))) ;
	}
	if (t<pow(2.0,-1.0) && t>=pow(3.0,-1.0)) {
		P[0] =-6.0*M_PI*R*sin(6.0*M_PI*(t-pow(3.0,-1.0))) ;
		P[1] = 6.0*M_PI*R*cos(6.0*M_PI*(t-pow(3.0,-1.0))) ;
	}
	if (t<=1.0 && t>=pow(2.0,-1.0)) {
		P[0] = 6.0*R*M_PI*sin(2.0*M_PI*(t-pow(2.0,-1.0))) ;
		P[1] =-6.0*R*M_PI*cos(2.0*M_PI*(t-pow(2.0,-1.0))) ;
	}
}

void KidneyRepSecDer(Point & P, double t) {
	//
	double R = 0.2 ;
	//
	if (t>=0 && t<pow(6.0,-1.0)) {
		P[0] =-pow(6.0*M_PI,2.0)*R*cos(6.0*M_PI*t) ;
		P[1] =-pow(6.0*M_PI,2.0)*R*sin(6.0*M_PI*t) ;
	}
	if (t<pow(3.0,-1.0) && t>=pow(6.0,-1.0)) {
		P[0] =-pow(6.0*M_PI,2.0)*R*cos(6.0*M_PI*(t-pow(6.0,-1.0))) ;
		P[1] =+pow(6.0*M_PI,2.0)*R*sin(6.0*M_PI*(t-pow(6.0,-1.0))) ;
	}
	if (t<pow(2.0,-1.0) && t>=pow(3.0,-1.0)) {
		P[0] =-pow(6.0*M_PI,2.0)*R*cos(6.0*M_PI*(t-pow(3.0,-1.0))) ;
		P[1] =-pow(6.0*M_PI,2.0)*R*sin(6.0*M_PI*(t-pow(3.0,-1.0))) ;
	}
	if (t<=1.0 && t>=pow(2.0,-1.0)) {
		P[0] = 3.0*pow(2.0*M_PI,2.0)*R*cos(2.0*M_PI*(t-pow(2.0,-1.0))) ;
		P[1] = 3.0*pow(2.0*M_PI,2.0)*R*sin(2.0*M_PI*(t-pow(2.0,-1.0))) ;
	}
}
	
#endif