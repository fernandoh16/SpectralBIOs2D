#ifndef FOURIERBOUNDARYREP_HEADERDEF 
#define FOURIERBOUNDARYREP_HEADERDEF

void R_Fourier(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] = cos(2.0 * M_PI * (double)(Index) * t) ;
		// P[1] = sin(2.0 * M_PI * (double)(Index) * t) ;
		P[0] = cos(2.0 * M_PI * t) ;
		P[1] = sin(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}

void R_der_Fourier(Point & P, double t, int Index) {
	if(Index == 0) {
		// P[0] =-sin(2.0 * M_PI * (double)(Index) * t) ;
		// P[1] = cos(2.0 * M_PI * (double)(Index) * t) ;
		P[0] =-2.0 * M_PI * sin(2.0 * M_PI * t) ;
		P[1] = 2.0 * M_PI * cos(2.0 * M_PI * t) ;
	}
	else {
		P[0] = 0.0 ;
		P[1] = 0.0 ;
	}
}
	
#endif