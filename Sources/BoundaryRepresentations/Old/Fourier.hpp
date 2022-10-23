#ifndef FOURIER_HEADERDEF 
#define FOURIER_HEADERDEF

void FourierPBR(Point & P, double t, int Index) {
	//
	double Value = 0.0 ;
	if(Index % 2 == 0) {
		Value = cos(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
	}
	else {
		Value = sin(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
	}
	P[0] = Value * cos(2.0 * M_PI * t) ;
	P[1] = Value * sin(2.0 * M_PI * t) ;
}

void FourierPBRDer(Point & P, double t, int Index) {
	double Value1 = 0.0 ;
	double Value2 = 0.0 ;
	if(Index % 2 == 0) {
		Value1 =                                       cos(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
		Value2 = -2.0 * M_PI * (double)(Index  )/2.0 * sin(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
	}
	else {
		Value1 =                                       sin(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
		Value2 = +2.0 * M_PI * (double)(Index+1)/2.0 * cos(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
	}
	P[0] = -2.0 * M_PI * Value1 * sin(2.0 * M_PI * t) + Value2 * cos(2.0 * M_PI * t) ;
	P[1] = +2.0 * M_PI * Value1 * cos(2.0 * M_PI * t) + Value2 * sin(2.0 * M_PI * t) ;
}

void FourierPBRSecDer(Point & P, double t, int Index) {
	double Value1 = 0.0 ;
	double Value2 = 0.0 ;
	double Value3 = 0.0 ;
	if(Index % 2 == 0) {
		Value1 =                                                cos(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
		Value2 =     -2.0 * M_PI * (double)(Index  )/2.0      * sin(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
		Value3 = -pow(2.0 * M_PI * (double)(Index  )/2.0,2.0) * cos(2.0 * M_PI * (double)(Index  )/2.0 * t) ;
	}
	else {
		Value1 =                                                sin(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
		Value2 =     +2.0 * M_PI * (double)(Index+1)/2.0      * cos(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
		Value3 = -pow(2.0 * M_PI * (double)(Index+1)/2.0,2.0) * sin(2.0 * M_PI * (double)(Index+1)/2.0 * t) ;
	}
	P[0] = Value3 * cos(2.0 * M_PI * t) - 4.0 * M_PI * Value2 * sin(2.0 * M_PI * t) ;
	P[1] = Value3 * sin(2.0 * M_PI * t) + 4.0 * M_PI * Value2 * cos(2.0 * M_PI * t) ;
	//
	P[0] = P[0] - pow(2.0 * M_PI,2.0) * Value1 * cos(2.0 * M_PI * t) ;
	P[1] = P[1] - pow(2.0 * M_PI,2.0) * Value1 * sin(2.0 * M_PI * t) ;
}

#endif