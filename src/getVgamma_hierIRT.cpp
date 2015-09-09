// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getVgamma_hierIRT(arma::cube &Vgamma,
					arma::mat &gammasigma,
					arma::mat &Ebb,
					const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &z,
                    const int NL,
                    const int NG
                    ) {
				
	signed int m, l;

#pragma omp parallel for private(m,l)
  	for(m=0; m < NG; m++){

		Vgamma.slice(m) = inv_sympd(gammasigma);

		for(l=0; l < NL; l++){
			if( g(i(l,0),0)==m ) Vgamma.slice(m) += Ebb(j(l,0),0) * trans(z.row(i(l,0))) * z.row(i(l,0));
		}

		// Note this yield C_inv, not C_m
		Vgamma.slice(m) = inv_sympd(Vgamma.slice(m));


  	}  // for(k=0; n < NJ; k++){

    return; 
} 
