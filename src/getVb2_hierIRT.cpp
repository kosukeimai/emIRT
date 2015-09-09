// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getVb2_hierIRT(arma::cube &Vb2,
					arma::mat &betasigma,
					arma::cube &Ex2x2,
					const arma::mat &i,
                    const arma::mat &j,
                    const int NL,
                    const int NJ
                    ) {
				
	// k loops over the NJ bills
	signed int k, l;

#pragma omp parallel for private(k,l)
  	for(k=0; k < NJ; k++){

		Vb2.slice(k) =  inv_sympd(betasigma);

		for(l=0; l < NL; l++){
			if( j(l,0)==k ) Vb2.slice(k) += Ex2x2.slice(i(l,0)) ;
		}

		// Note this yield B_inv, not B
		Vb2.slice(k) = inv_sympd(Vb2.slice(k));

  	}  // for(k=0; n < NJ; k++){

    return; 
} 
