// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEgamma_hierIRT(arma::mat &Egamma,
					const arma::cube &Vgamma,
					const arma::mat &gammasigma,
					const arma::mat &gammamu,
					const arma::mat &g,
					const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &z,
                    const arma::mat &Eb,
                    const arma::mat &Ebb,
                    const arma::mat &Eystar,
                    const arma::mat &Ea,
                    const arma::mat &Eta,
                    const int NL,
                    const int NG
                    ) {

	signed int m, l;

#pragma omp parallel for private(m,l)
  	for(m=0; m < NG; m++){

		Egamma.row(m) = trans(inv_sympd(gammasigma) * gammamu);

		for(l=0; l < NL; l++){
			if( g(i(l,0),0)==m ){

				Egamma.row(m) += z.row(i(l,0)) * (  Eb(j(l,0),0)*(Eystar(l,0) - Ea(j(l,0),0)) - Ebb(j(l,0),0)*Eta(i(l,0),0) )  ;

			}
		}

		// Note this yield C_inv, not C_m
		Egamma.row(m) = trans( Vgamma.slice(m) * trans(Egamma.row(m)) );

  	}  // for(k=0; n < NJ; k++){

    return; 
} 

