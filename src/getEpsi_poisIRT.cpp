// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEpsi(arma::mat &Epsi,
				   const arma::mat &Exfull,
				   const arma::mat &onesNK,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &i,
                   const arma::mat &Ealpha,
                   const arma::mat &Ebeta,
                   const double psi_mu,
                   const double psi_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int j;
	double denom, numer;

#pragma omp parallel for private(denom, numer)
	for(j=0; j < NJ; j++){

		denom = as_scalar( exi.row(j) * onesNK )  + 1/psi_sigma;
		numer = as_scalar(y.row(j)*onesNK - denom + 1/psi_sigma - exi.row(j)*Ealpha + exixi.row(j)*onesNK - exi.row(j)*Exfull*Ebeta(j,0)) + psi_mu/psi_sigma;
	
		Epsi(j,0) = numer/denom;
		
	}

    return;

}

