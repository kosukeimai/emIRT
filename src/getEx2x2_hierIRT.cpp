// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEx2x2_hierIRT(arma::cube &curEx2x2,
					const arma::mat &curEx2,
                    const arma::cube &Egg,
                    const arma::mat &Eg,
                    const arma::mat &Eta,
                    const arma::mat &Eta2,
                    const arma::mat &g,
                    const arma::mat &z,
                    const int NI
                    ) {


	signed int n;	//indexed each i

    // Main Calculation
#pragma omp parallel for
  	for(n=0; n < NI; n++){

		curEx2x2(0,1,n) = curEx2(n,1);
		curEx2x2(1,0,n) = curEx2(n,1);

		// Independent approximation
		//curEx2x2(1,1,n) = curEx2(n,1)*curEx2(n,1);

	 	curEx2x2(1,1,n) = accu( z.row(n) * Egg.slice(g(n,0)) * trans(z.row(n)) ) + 2*accu(Eg.row(g(n,0))*trans(z.row(n))*Eta.row(n)) + Eta2(n,0) ;


  	}  // for(n=0; n < NI; n++){

    return; 
} 

