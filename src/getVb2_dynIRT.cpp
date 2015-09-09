// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getVb2_dynIRT(arma::cube &Vb2,
				 const arma::cube &Ex2x2,
                 const arma::mat &sigma,
                 const int T
                 ) {

	int t;

#pragma omp parallel for
	for(t=0; t<T; t++){

		Vb2.slice(t) = inv_sympd(inv_sympd(sigma) + Ex2x2.slice(t)) ;

	}

    return;

}
