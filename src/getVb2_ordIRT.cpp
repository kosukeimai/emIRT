// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::cube getVb2_ordIRT(const arma::mat &Ex,
                 const arma::mat &Ex2x2,
                 const arma::mat &sigma,
                 const arma::mat &Edd,
                 const int nJ
                 ) {

	int j;
	arma::cube Vb2(2, 2, nJ);

#pragma omp parallel for
 	for(j=0; j < nJ; j++){

	    Vb2.slice(j) = inv_sympd(inv_sympd(sigma) + Edd(j,0)*Ex2x2 );

	}

    return(Vb2) ;
}
