// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEalpha(arma::mat &Ealpha,
				   arma::mat &Ex,
				   const arma::mat &onesNJ,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &i,
                   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const double alpha_mu,
                   const double alpha_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int j, k;
	double denom_alpha, numer_alpha;
	arma::mat Ebb(NJ,1);
	arma::mat Ebp(NJ,1);

#pragma omp parallel for 
	for(j=0; j < NJ; j++){
		Ebb(j,0) = Ebeta(j,0)*Ebeta(j,0);
		Ebp(j,0) = Ebeta(j,0)*Epsi(j,0);
	}

#pragma omp parallel for private(denom_alpha, numer_alpha)
	for(k=0; k < NK; k++){

		denom_alpha = as_scalar( trans(exi.col(k)) * onesNJ )  + 1/alpha_sigma;
		numer_alpha = as_scalar(-1*denom_alpha + trans(y.col(k))*onesNJ - trans(exi.col(k))*Epsi + trans(xi.col(k))*exi.col(k) - trans(exi.col(k))*Ebeta*Ex(i(k,0),0) ) + 1/alpha_sigma + alpha_mu/alpha_sigma;
		Ealpha(k,0) = numer_alpha/denom_alpha;

	}
		
    return;

}

