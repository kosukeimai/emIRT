// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEb2_hierIRT(arma::mat &Eb2,
					const arma::cube &Vb2,
					const arma::mat &betasigma,
					const arma::mat &betamu,
					const arma::mat &Eystar,
					const arma::mat &Ex2,
					const arma::mat &i,
                    const arma::mat &j,
                    const int NL,
                    const int NJ
                    ) {

	// k loops over the NJ bills
	signed int k, l;

#pragma omp parallel for private(k,l)
  	for(k=0; k < NJ; k++){

		Eb2.row(k) =  trans(inv_sympd(betasigma) * betamu);

		for(l=0; l < NL; l++){
			if( j(l,0)==k ) Eb2.row(k) += Eystar(l,0) * Ex2.row(i(l,0));
		}

		// Multiplying by variance and retransposing 
		Eb2.row(k) = trans(Vb2.slice(k) * trans(Eb2.row(k)));

  	}  // for(k=0; n < NJ; k++){

    return; 
} 
