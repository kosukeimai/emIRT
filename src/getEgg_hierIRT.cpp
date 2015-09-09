// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEgg_hierIRT(arma::cube &Egg,
					const arma::mat &Egamma,
					const arma::cube &Vgamma,
                    const int NG
                    ) {

	signed int m;

#pragma omp parallel for
  	for(m=0; m < NG; m++){

		Egg.slice(m) = Vgamma.slice(m) + trans(Egamma.row(m)) * Egamma.row(m);

  	}  // for(k=0; n < NJ; k++){

    return; 
} 
