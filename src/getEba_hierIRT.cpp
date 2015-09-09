// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEba_hierIRT(arma::mat &Eba,
					const arma::mat &Eb,
					const arma::mat &Ea,
					const arma::cube &Vb2,
                    const int NJ
                    ) {

	signed int k;

#pragma omp parallel for
  	for(k=0; k < NJ; k++){

		Eba(k,0) = Eb(k,0)*Ea(k,0) + Vb2(0,1,k);

  	}  // for(k=0; n < NJ; k++){

    return; 
} 
