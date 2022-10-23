#ifndef BOOMERANG_HEADERDEF 
#define BOOMERANG_HEADERDEF

void BoomerangRep(Point & P, double t) {
	P[0] = 0.5 + cos(2.0*M_PI*t) + 0.65 * cos(4.0*M_PI*t) - 0.65 ;
	P[1] =     + 1.5*sin(2.0*M_PI*t) ;
	//
	double R = 0.25 ;
	P = P * R ;
}

void BoomerangRepDer(Point & P, double t) {
	P[0] = -2.0*M_PI*sin(2.0*M_PI*t)-0.65*4.0*M_PI*sin(4.0*M_PI*t) ;
	P[1] =  3.0*M_PI*cos(2.0*M_PI*t) ;	
	//
	double R = 0.25 ;
	P = P * R ;
}

void BoomerangRepSecDer(Point & P, double t) {
	P[0] = -pow(2.0*M_PI,2.0)*cos(2.0*M_PI*t) - pow(4.0*M_PI,2.0)*0.65*cos(4.0 * M_PI * t) ;
	P[1] = -1.5*pow(2.0*M_PI,2.0)*sin(2.0*M_PI*t) ;
	//
	double R = 0.25 ;
	P = P * R ;
}

void BoomerangRep2(Point & P, double t) {
	P[0] = 5.0 * cos(2.0*M_PI*t) - 3.25 * cos(4.0*M_PI*t) ;
	P[1] = 7.5 * sin(2.0*M_PI*t) ;
	//
	/*
	double R = 1.0 ;
	P = P * R ;
	*/
}

void BoomerangRep2Der(Point & P, double t) {
	P[0] = -10.0 * M_PI * sin(2.0*M_PI*t) + 3.25*4.0*M_PI*sin(4.0*M_PI*t) ;
	P[1] =  15.0 * M_PI * cos(2.0*M_PI*t) ;	
	//
	/*
	double R = 1.0 ;
	P = P * R ;
	*/
}

void BoomerangRep2SecDer(Point & P, double t) {
	P[0] = -5.0 * pow(2.0*M_PI,2.0)*cos(2.0*M_PI*t) + pow(4.0*M_PI,2.0) * 3.25 * cos(4.0 * M_PI * t) ;
	P[1] = -7.5 * pow(2.0*M_PI,2.0)*sin(2.0*M_PI*t) ;
	//
	/*
	double R = 0.4 ;
	P = P * R ;
	*/
}

#endif