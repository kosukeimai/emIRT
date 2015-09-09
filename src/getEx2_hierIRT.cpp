// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEx2_hierIRT(arma::mat &curEx2,
					const arma::mat &Egamma,
                    const arma::mat &Eta,
                    const arma::mat &g,
                    const arma::mat &z,
                    const int NI
                    ) {

	signed int n;	//indexed each i

    // Main Calculation
#pragma omp parallel for
  	for(n=0; n < NI; n++){

		curEx2(n,1) = accu(Egamma.row(g(n,0)) % z.row(n)) + Eta(n,0);

  	}  // for(n=0; n < NI; n++){

    return; 
} 
