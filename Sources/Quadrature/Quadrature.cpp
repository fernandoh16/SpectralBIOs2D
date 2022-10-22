#include <cmath>
#include <string>
#include <iostream>
#include <complex>
#include <cassert>

#include "Quadrature.hpp"

Quadrature::Quadrature() {

}
void Quadrature::SetQuadrature(string Dimension_, int  NumberofPoints_, int NumberofCycles_) {
	assert((Dimension_.compare("1D")==0) || (Dimension_.compare("2D")==0));
	Dimension = Dimension_; 
	NumberofPoints = NumberofPoints_;
	NumberofCycles = NumberofCycles_ ;
	if (Dimension.compare("1D")==0) {
		QuadraturePoints.reserve(NumberofPoints)  ;
		QuadratureWeights.reserve(NumberofPoints) ;
		if (NumberofPoints==2) {
			//
			QuadraturePoints[0] = -0.57735;
			QuadraturePoints[1] = +0.57735;
			//
			QuadratureWeights[0] = +1.00000;
			QuadratureWeights[1] = +1.00000;
		}
	    if (NumberofPoints==3) {
			//
			QuadraturePoints[0] = -0.77459;
			QuadraturePoints[1] = +0.00000;
			QuadraturePoints[2] = +0.77459;
			//
			QuadratureWeights[0] = +0.55555;
			QuadratureWeights[1] = +0.88888;
			QuadratureWeights[2] = +0.55555;
		}
		if (NumberofPoints==4) {
			//
			QuadraturePoints[0] = -0.8611363116;
			QuadraturePoints[1] = -0.3399810436;
			QuadraturePoints[2] = +0.3399810436;
			QuadraturePoints[3] = +0.8611363116;
			//
			QuadratureWeights[0] = +0.3478548451;
			QuadratureWeights[1] = +0.6521451549;
			QuadratureWeights[2] = +0.6521451549;
			QuadratureWeights[3] = +0.3478548451;
		}
		if (NumberofPoints==12) {
			//
			QuadraturePoints[0]  = -0.981560634246719;
			QuadraturePoints[1]  = -0.904117256370475;
			QuadraturePoints[2]  = -0.769902674194305;
			QuadraturePoints[3]  = -0.587317954286617;
			QuadraturePoints[4]  = -0.367831498998180;
			QuadraturePoints[5]  = -0.125233408511469;
			QuadraturePoints[6]	= +0.125233408511469;
			QuadraturePoints[7]  = +0.367831498998180;
			QuadraturePoints[8]  = +0.587317954286617;
			QuadraturePoints[9]  = +0.769902674194305;
			QuadraturePoints[10] = +0.904117256370475;
			QuadraturePoints[11] = +0.981560634246719;
			//
			QuadratureWeights[0]  = +0.047175336386512;
			QuadratureWeights[1]  = +0.106939325995318;
			QuadratureWeights[2]  = +0.160078328543346;
			QuadratureWeights[3]  = +0.203167426723066;
			QuadratureWeights[4]  = +0.233492536538355;
			QuadratureWeights[5]  = +0.249147045813403;
			QuadratureWeights[6]  = +0.249147045813403;
			QuadratureWeights[7]  = +0.233492536538355;
			QuadratureWeights[8]  = +0.203167426723066;
			QuadratureWeights[9]  = +0.160078328543346;
			QuadratureWeights[10] = +0.106939325995318;
			QuadratureWeights[11] = +0.047175336386512;
		}
		if (NumberofPoints==24) {
			//
			QuadraturePoints[0]  = -0.995187219997021360180;
			QuadraturePoints[1]  = -0.974728555971309498198;
			QuadraturePoints[2]  = -0.938274552002732758524;
			QuadraturePoints[3]  = -0.886415527004401034213;
			QuadraturePoints[4]  = -0.820001985973902921954;
			QuadraturePoints[5]  = -0.740124191578554364244;
			QuadraturePoints[6]	 = -0.648093651936975569252;
			QuadraturePoints[7]  = -0.545421471388839535658;
			QuadraturePoints[8]  = -0.433793507626045138487;
			QuadraturePoints[9]  = -0.315042679696163374387;
			QuadraturePoints[10] = -0.191118867473616309159;
			QuadraturePoints[11] = -0.064056892862605626085;
			QuadraturePoints[12] = +0.064056892862605626085;
			QuadraturePoints[13] = +0.191118867473616309159;
			QuadraturePoints[14] = +0.315042679696163374387;
			QuadraturePoints[15] = +0.433793507626045138487;
			QuadraturePoints[16] = +0.545421471388839535658;
			QuadraturePoints[17] = +0.648093651936975569252;
			QuadraturePoints[18] = +0.740124191578554364244;
			QuadraturePoints[19] = +0.820001985973902921954;
			QuadraturePoints[20] = +0.886415527004401034213;
			QuadraturePoints[21] = +0.938274552002732758524;
			QuadraturePoints[22] = +0.974728555971309498198;
			QuadraturePoints[23] = +0.995187219997021360180;
			//
			QuadratureWeights[0]  = +0.012341229799987199547;
			QuadratureWeights[1]  = +0.028531388628933663181;
			QuadratureWeights[2]  = +0.044277438817419806169;
			QuadratureWeights[3]  = +0.059298584915436780746;
			QuadratureWeights[4]  = +0.073346481411080305734;
			QuadratureWeights[5]  = +0.086190161531953275917;
			QuadratureWeights[6]  = +0.097618652104113888270;
			QuadratureWeights[7]  = +0.107444270115965634783;
			QuadratureWeights[8]  = +0.115505668053725601353;
			QuadratureWeights[9]  = +0.121670472927803391204;
			QuadratureWeights[10] = +0.125837456346828296121;
			QuadratureWeights[11] = +0.127938195346752156974;
			QuadratureWeights[12] = +0.127938195346752156974;
			QuadratureWeights[13] = +0.125837456346828296121;
			QuadratureWeights[14] = +0.121670472927803391204;
			QuadratureWeights[15] = +0.115505668053725601353;
			QuadratureWeights[16] = +0.107444270115965634783;
			QuadratureWeights[17] = +0.097618652104113888270;
			QuadratureWeights[18] = +0.086190161531953275917;
			QuadratureWeights[19] = +0.073346481411080305734;
			QuadratureWeights[20] = +0.059298584915436780746;
			QuadratureWeights[21] = +0.044277438817419806169;
			QuadratureWeights[22] = +0.028531388628933663181;
			QuadratureWeights[23] = +0.012341229799987199547;
		}
	}
	if (Dimension.compare("2D")==0) {
		QuadraturePoints.reserve(2*NumberofPoints)  ;
		QuadratureWeights.reserve(NumberofPoints) ;
		if (NumberofPoints==4) {
			//
			double aux = sqrt((double)(1)/(double)(3));	
			QuadraturePoints[0] = -aux; QuadraturePoints[1] = -aux;
			QuadraturePoints[2] = +aux; QuadraturePoints[3] = -aux;
			QuadraturePoints[4] = -aux; QuadraturePoints[5] = +aux;
			QuadraturePoints[6] = +aux; QuadraturePoints[7] = +aux;
			//
			QuadratureWeights[0] = +0.25000;
			QuadratureWeights[1] = +0.25000;
			QuadratureWeights[2] = +0.25000;
			QuadratureWeights[3] = +0.25000;	
		}
		if (NumberofPoints==9) {	
			//
			double aux = sqrt((double)(3)/(double)(5));		
			QuadraturePoints[0]  = -aux; QuadraturePoints[1]  = -aux;
			QuadraturePoints[2]  = -aux; QuadraturePoints[3]  = +aux;
			QuadraturePoints[4]  = +aux; QuadraturePoints[5]  = -aux;
			QuadraturePoints[6]  = +aux; QuadraturePoints[7]  = +aux;
			QuadraturePoints[8]  = -aux; QuadraturePoints[9]  = +0.0;
			QuadraturePoints[10] = +aux; QuadraturePoints[11] = +0.0;
			QuadraturePoints[12] = +0.0; QuadraturePoints[13] = -aux;
			QuadraturePoints[14] = +0.0; QuadraturePoints[15] = +aux;
			QuadraturePoints[16] = +0.0; QuadraturePoints[17] = +0.0;
			//
			QuadratureWeights[0] = +(double)(25)/(double)(324);
			QuadratureWeights[1] = +(double)(25)/(double)(324);
			QuadratureWeights[2] = +(double)(25)/(double)(324);
			QuadratureWeights[3] = +(double)(25)/(double)(324);
			QuadratureWeights[4] = +(double)(10)/(double)(81);
			QuadratureWeights[5] = +(double)(10)/(double)(81);
			QuadratureWeights[6] = +(double)(10)/(double)(81);
			QuadratureWeights[7] = +(double)(10)/(double)(81);
			QuadratureWeights[8] = +(double)(16)/(double)(81);	
		}
	}
}

double Quadrature::ComputeQuadrature(double (*p_func)(double x, void * Input), void * Input_, double x1, double x2) {
	double QuadratureValue   ;
	double QuadraturePoint   ;
	double QuadratureWeight  ;
	double ReturnValue = 0.0 ;
	double xA ;
	double xB ;
	
	for(int Cycle = 0; Cycle<NumberofCycles; Cycle ++) {
		QuadratureValue = 0.0 ;
		xA = x1 + (double)(Cycle  )/(double)(NumberofCycles) * (x2-x1) ;
		xB = x1 + (double)(Cycle+1)/(double)(NumberofCycles) * (x2-x1) ;
		for(int i = 0; i<NumberofPoints; i++) {
			QuadraturePoint  = (QuadraturePoints[i])*(xB-xA)/2.0 + (xB+xA)/2.0  ;
			QuadratureWeight = QuadratureWeights[i] ;
			QuadratureValue  = QuadratureValue + QuadratureWeight * (*p_func)(QuadraturePoint,Input_) ;
		}
		QuadratureValue = QuadratureValue*(xB-xA)/2.0 ;
		ReturnValue = ReturnValue + QuadratureValue ;
	}
	return ReturnValue ;
}

double Quadrature::ComputeQuadrature(double (*p_func)(double x, double y, void * Input), void * Input_, double x1, double x2, double y1, double y2) {
	double QuadratureValue = 0.0 ;
	double QuadraturePoint1 ;
	double QuadraturePoint2 ;
	double QuadratureWeight ;
	for(int i = 0; i<NumberofPoints; i++) {
		QuadraturePoint1  = (QuadraturePoints[2*i    ])*(x2-x1)/2.0 + (x2+x1)/2.0 ;
		QuadraturePoint2  = (QuadraturePoints[2*i + 1])*(y2-y1)/2.0 + (y2+y1)/2.0 ;
		QuadratureWeight  = QuadratureWeights[i] ;
		QuadratureValue   = QuadratureValue + QuadratureWeight*(*p_func)(QuadraturePoint1,QuadraturePoint2,Input_) ;
	}
	QuadratureValue = QuadratureValue*(x2-x1)*(y2-y1) ;
	return QuadratureValue ;
}

complex < double > Quadrature::ComputeQuadrature(complex < double > (*p_func)(double x, double y, void * Input), void * Input_, double x1, double x2, double y1, double y2) {
	complex < double > QuadratureValue = complex < double > (0.0,0.0) ;
	double QuadraturePoint1 ;
	double QuadraturePoint2 ;
	double QuadratureWeight ;
	complex < double > ReturnValue = complex < double > (0.0,0.0) ;
	double xA1 ; double xA2 ;
	double yB1 ; double yB2 ;
	for(int Cycle1 = 0; Cycle1<NumberofCycles; Cycle1++) {  
		for(int Cycle2 = 0; Cycle2<NumberofCycles; Cycle2++) {
			QuadratureValue = complex < double > (0.0,0.0) ;
			xA1 = x1 + (double)(Cycle1  )/(double)(NumberofCycles) * (x2-x1) ;
			xA2 = x1 + (double)(Cycle1+1)/(double)(NumberofCycles) * (x2-x1) ;
			yB1 = y1 + (double)(Cycle2  )/(double)(NumberofCycles) * (y2-y1) ;
			yB2 = y1 + (double)(Cycle2+1)/(double)(NumberofCycles) * (y2-y1) ;
			for(int i = 0; i<NumberofPoints; i++) {
				QuadraturePoint1  = (QuadraturePoints[2*i    ])*(xA2-xA1)/2.0 + (xA2+xA1)/2.0 ;
				QuadraturePoint2  = (QuadraturePoints[2*i + 1])*(yB2-yB1)/2.0 + (yB2+yB1)/2.0 ;
				QuadratureWeight  = QuadratureWeights[i] ;
				QuadratureValue   = QuadratureValue + QuadratureWeight*(*p_func)(QuadraturePoint1,QuadraturePoint2,Input_) ;
			}
			QuadratureValue = QuadratureValue*(xA2-xA1)*(yB2-yB1) ;
			ReturnValue = ReturnValue + QuadratureValue ;
		}
	}
	return ReturnValue ;
}

int Quadrature::GetNumberOfPoints() {
	return NumberofPoints ;
}

double Quadrature::GetQuadraturePoint(int i, double x1, double x2) {
	return (QuadraturePoints[i])*(x2-x1)/2.0 + (x2+x1)/2.0 ;
}

double Quadrature::GetQuadratureWeight(int i, double x1, double x2) {
	return (QuadratureWeights[i])*(x2-x1)/2.0 ;
}