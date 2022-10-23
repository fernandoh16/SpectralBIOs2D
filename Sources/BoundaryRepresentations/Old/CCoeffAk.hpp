#ifndef COEFFAK_HEADERDEF 
#define COEFFAK_HEADERDEF

// Generate Coefficients of the Form j^{—zeta}
void ConstructionCoeffAk(vector<double> & CoeffAK, double Zeta) {
	for(int i = 1; i<=CoeffAK.size(); i++) {
		CoeffAK[i-1] = pow((double)(i),-Zeta) ;
	}
	for(int i = 1; i<=CoeffAK.size(); i++) {
		if(i%4 == 0) {
			CoeffAK[i-1] = pow((double)(i)/(double)(4),-Zeta) ;
		}
		if((i+1)%4 == 0) {
			CoeffAK[i-1] = pow((double)(i+1)/(double)(4),-Zeta) ;
		}
		if((i+2)%4 == 0) {
			CoeffAK[i-1] = pow((double)(i+2)/(double)(4),-Zeta) ;
		}
		else if((i+3)%4 == 0) {
			CoeffAK[i-1] = pow((double)(i+3)/(double)(4),-Zeta) ;
		}
	}
	//
	/*
	div_t  div_res ;
	//
	double Index_ ;
	for(int i = 1; i<=CoeffAK.size(); i++) {
		div_res = div(i+4,4) ;  
		Index_ = (double)(div_res.quot) ;
		if((i+3)%4 == 0) {
			CoeffAK[i-1] = pow(Index_,-Zeta) ;
		}
		else if((i+2)%4 == 0) {
			CoeffAK[i-1] = pow(Index_,-Zeta) ;
		}
		else if((i+1)%4 == 0) {
			CoeffAK[i-1] = pow(Index_,-Zeta) ;
		}
		else if((i+0)%4 == 0) {
			CoeffAK[i-1] = pow(Index_,-Zeta) ;
		}
		cout << Index_ << "\n" ;
	}
	*/
}
	
#endif