// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVeta_hierIRT(const arma::mat &i,
                    const arma::mat &j,
                    const arma::mat &g,
                    const arma::mat &curEsigma,
                    const arma::mat &curEbb,
                    const int NI,
                    const int NL
                    ) {

	signed int n, l;
	arma::mat Veta(NI,1);
	double tmp;

    // Main Calculation
#pragma omp parallel for private(n,l,tmp)
  	for(n=0; n < NI; n++){

		tmp = 1/curEsigma(g(n,0),0);

		for(l=0; l < NL; l++){
			if( i(l,0)==n ) tmp += curEbb(j(l,0),0);
		}

		//This gives A_inv, instead of A.  A_inv is actually the variance
		Veta(n,0) = 1/tmp;


  	}  // for(n=0; n < NI; n++){

    return(Veta); 
} 
