#ifndef BAYESIANINVERSION_HPP
#define BAYESIANINVERSION_HPP

#include <direct_product.hpp>

template <class F, class P, class M, class Q, class S, class V> 
class BayesianInversion {
private:
    F & ForwardModel ;
    P & Potential    ;
    M & Measurements ;
    Q & QuantityOfInterest ;
public:
    BayesianInversion(F & ForwardModel_, P & Potential_, M & Measurements_, Q & QuantityOfInterest_) ;
    vectorspace::direct_product<V, V, double> operator()(S y) ;
} ;

template <class F, class P, class M, class Q, class S, class V> 
BayesianInversion<F, P, M, Q, S, V>::BayesianInversion(F & ForwardModel_, P & Potential_, M & Measurements_, Q & QuantityOfInterest_) : 
	ForwardModel(ForwardModel_), 
	Potential(Potential_), 
	Measurements(Measurements_),
	QuantityOfInterest(QuantityOfInterest_) {
}

template <class F, class P, class M, class Q, class S, class V> 
vectorspace::direct_product<V, V, double> BayesianInversion<F, P, M, Q, S, V>::operator()(S y) {
	auto 		 QoI = QuantityOfInterest(y) ;
	auto Observation = ForwardModel(y) ;
    double Theta = exp(-Potential(Observation,Measurements)) ;
    return vectorspace::direct_product<V, V, double>(QoI, Theta * QoI, Theta) ;
}

#endif

				