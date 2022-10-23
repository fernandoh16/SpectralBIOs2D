#ifndef CIRCLE_HEADERDEF 
#define CIRCLE_HEADERDEF

void CircleRep(Point & P, double t) {
	P[0] = cos(2.0 * M_PI * t) ;
	P[1] = sin(2.0 * M_PI * t) ;
}

void CircleRepDer(Point & P, double t) {
	P[0] = -2.0 * M_PI * sin(2.0 * M_PI * t) ;
	P[1] = +2.0 * M_PI * cos(2.0 * M_PI * t) ;
}

void CircleRepSecDer(Point & P, double t) {
	P[0] = -pow(2.0 * M_PI,2.0) * cos(2.0 * M_PI * t) ;
	P[1] = -pow(2.0 * M_PI,2.0) * sin(2.0 * M_PI * t) ;
}

void CircleRep2(Point & P, double t) {
	double R = 0.3 ;
	P[0] = R * cos(2.0 * M_PI * t) ;
	P[1] = R * sin(2.0 * M_PI * t) ;
}

void CircleRep2Der(Point & P, double t) {
	double R = 0.3 ;
	P[0] = -2.0 * M_PI * R * sin(2.0 * M_PI * t) ;
	P[1] = +2.0 * M_PI * R * cos(2.0 * M_PI * t) ;
}

void CircleRep2SecDer(Point & P, double t) {
	double R = 0.3 ;
	P[0] = -pow(2.0 * M_PI,2.0) * R * cos(2.0 * M_PI * t) ;
	P[1] = -pow(2.0 * M_PI,2.0) * R * sin(2.0 * M_PI * t) ;
}

#endif