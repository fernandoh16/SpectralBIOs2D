#ifndef KIDNEY_HEADERDEF 
#define KIDNEY_HEADERDEF
/*
#include <stdlib.h> 
#include <cmath> 

using namespace std ;
*/
void Kidney_rep(Point & P, double t, int Index) {
	//
	div_t  div_res ;
	div_res = div(Index+4,4) ;  
	//
	double Index_ = (double)(div_res.quot) ;
	//
	double   eta = -4.0 ;
	double sigma =  0.15 ;
	double R = 0.04/0.15 ;
	//
	if(Index == 0) {
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
	else if((Index+3) % 4 == 0) {
		P[0] = sigma * pow(Index_,eta) * cos(2.0 * M_PI * Index_ * t) ;
		P[1] = 0.0 ;
	}
	else if((Index+2) % 4 == 0) {
		P[0] = sigma * pow(Index_,eta) * sin(2.0 * M_PI * Index_ * t) ;
		P[1] = 0.0 ;
	}
	else if((Index+1) % 4 == 0) {
		P[0] = 0.0 ;
		P[1] = sigma * pow(Index_,eta) * cos(2.0 * M_PI * Index_ * t) ;
	}
	else if((Index+0) % 4 == 0) {
		P[0] = 0.0 ;
		P[1] = sigma * pow(Index_,eta) * sin(2.0 * M_PI * Index_ * t) ;
	}
}

void Kidney_rep_der(Point & P, double t, int Index) {
	//
	div_t  div_res ;
	div_res = div(Index+4,4) ;  
	//
	double Index_ = (double)(div_res.quot) ;
	//
	double   eta = -4.0 ;
	double sigma =  0.125 ;
	double R = 0.2 ;
	//
	if(Index == 0) {
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
			P[0] = 6.0*R*sin(2.0*M_PI*(t-pow(2.0,-1.0))) ;
			P[1] =-6.0*R*cos(2.0*M_PI*(t-pow(2.0,-1.0))) ;
		}
		
	}
	else if((Index+3) % 4 == 0) {
		P[0] = - 2.0 * M_PI * sigma * pow(Index_,eta) * Index_ * sin(2.0 * M_PI * Index_ * t) ;
		P[1] = 0.0 ;
	}
	else if((Index+2) % 4 == 0) {
		P[0] = + 2.0 * M_PI * sigma * pow(Index_,eta) * Index_ * cos(2.0 * M_PI * Index_ * t) ;
		P[1] = 0.0 ;
	}
	else if((Index+1) % 4 == 0) {
		P[0] = 0.0 ;
		P[1] = - 2.0 * M_PI * sigma * pow(Index_,eta) * Index_ * sin(2.0 * M_PI * Index_ * t) ;
	}
	else if((Index+0) % 4 == 0) {
		P[0] = 0.0 ;
		P[1] = + 2.0 * M_PI * sigma * pow(Index_,eta) * Index_ * cos(2.0 * M_PI * Index_ * t) ;
	}
}
	
#endif