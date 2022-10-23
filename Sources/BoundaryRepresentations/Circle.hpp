#ifndef FOURIERBOUNDARYREP_HEADERDEF 
#define FOURIERBOUNDARYREP_HEADERDEF

void R_Fourier(Point & P, double t, int Index) {
	if(Index == 0) {
		P[0] = 0.05 * cos(2.0 * M_PI * t) ;
		P[1] = 0.5 * sin(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}

void R_der_Fourier(Point & P, double t, int Index) {
	if(Index == 0) {
		P[0] =-0.05 * 2.0 * M_PI * sin(2.0 * M_PI * t) ;
		P[1] = 0.5 * 2.0 * M_PI * cos(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}

void R_sec_der_Fourier(Point & P, double t, int Index) {
	if(Index == 0) {
		P[0] =-0.05 * pow(2.0 * M_PI,2.0) * cos(2.0 * M_PI * t) ;
		P[1] =-0.5 * pow(2.0 * M_PI,2.0) * sin(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}


#endif