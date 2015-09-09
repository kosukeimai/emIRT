// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEsigma_hierIRT(arma::mat &Esigma,
					const arma::mat &Eta2,
					const arma::mat &sigmav,
					const arma::mat &sigmas,
					const arma::mat &g,
                    const int NG,
                    const int NI
                    ) {

	signed int m, n;
	double numerator, denominator;

#pragma omp parallel for private(n,m,denominator,numerator)
  	for(m=0; m < NG; m++){

		denominator = sigmav(0,0);
		numerator = sigmas(0,0);

		for(n=0; n < NI; n++){

			if( g(n,0)==m ){
				denominator += 1;
				numerator += Eta2(n,0);
			}

		}

		// Note this yields E[sigma^2], not E[sigma^{-2}]
		Esigma(m,0) = numerator/denominator;

  	}  // for(m=0; m < NG; m++){

    return; 
} 

