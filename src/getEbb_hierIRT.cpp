// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEbb_hierIRT(arma::mat &Ebb,
					const arma::mat &Eb,
					const arma::cube &Vb2,
                    const int NJ
                    ) {

	signed int k;

#pragma omp parallel for
  	for(k=0; k < NJ; k++){

		Ebb(k,0) = Eb(k,0)*Eb(k,0) + Vb2(1,1,k);

  	}  // for(k=0; n < NJ; k++){

    return; 
} 
