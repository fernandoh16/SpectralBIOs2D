#ifndef BOOMERANG_HEADERDEF 
#define BOOMERANG_HEADERDEF
/*
#include <stdlib.h> 
#include <cmath> 

using namespace std ;
*/
void Boomerang(Point & P, double t, int Index) {
	//
	div_t  div_res ;
	div_res = div(Index+4,4) ;  
	//
	double Index_ = (double)(div_res.quot) ;
	//
	double   eta = -4.0 ;
	double sigma =  0.2 ;
	//
	double F = 1.0 ;
	//
	if(Index == 0) {
		/*
		P[0] = 0.2*cos(2.0 * M_PI * t) + 0.5 * cos(4.0 * M_PI * t) - 0.65 ;
		P[1] = 1.5 *0.5/0.65* sin(2.0 * M_PI * t) ;
		P = P * F ;
		*/
		P[0] = F * cos(2.0 * M_PI * t) ;
		P[1] = F * sin(2.0 * M_PI * t) ;
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

void Boomerang_der(Point & P, double t, int Index) {
	//
	div_t  div_res ;
	div_res = div(Index+4,4) ;  
	//
	double Index_ = (double)(div_res.quot) ;
	//
	double   eta = -4.0 ;
	double sigma =  0.2 ;
	//
	double F = 1.0 ;
	//
	if(Index == 0) {	
		/*
		P[0] = -2.0 * M_PI * sin(2.0 * M_PI * t) - 4.0 * M_PI * 1.5 * sin(4.0 * M_PI * t) ;
		P[1] = 1.5 * 2.0 * M_PI * cos(2.0 * M_PI * t) ;
		P = P * F ;
		*/
		P[0] = -2.0 * M_PI * F * sin(2.0 * M_PI * t) ;
		P[1] = +2.0 * M_PI * F * cos(2.0 * M_PI * t) ;
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
/*
void Boomerang_sec_der(Point & P, double t, int Index) {
	//
	div_t  div_res ;
	div_res = div(Index+4,4) ;  
	//
	double Index_ = (double)(div_res.quot) ;
	//
	double   eta = -2.0 ;
	double sigma =  0.2 ;
	//
	if(Index == 0) {
		P[0] = -pow(2.0 * M_PI,2.0) * cos(2.0 * M_PI * t) - pow(4.0 * M_PI,2.0) * 0.65 * cos(4.0 * M_PI * t) ;
		P[1] = -1.5 * pow(2.0 * M_PI,2.0) * sin(2.0 * M_PI * t) ;
	}
	else if((Index+3) % 4 == 0) {
		P[0] = - 2.0 * M_PI * sigma * pow(Index_,eta+2) * sin(2.0 * M_PI * Index_ * t) ;
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
*/
	
#endif