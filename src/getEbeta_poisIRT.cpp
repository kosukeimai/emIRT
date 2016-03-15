// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEbeta(arma::mat &Ebeta,
				   const arma::mat &Vbeta,
				   const arma::mat &Ealpha,
				   const arma::mat &onesNK,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &Epsi,
                   const arma::mat &Exfull,
                   const arma::mat &Vx,
                   const arma::mat &i,
                   const double beta_mu,
                   const double beta_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int j, k;
	double numer;

	arma::mat Exa(NK,1);

#pragma omp parallel for
	for(k=0; k<NK; k++){
		Exa(k,0) = Exfull(k,0)*Ealpha(k,0);
	}

#pragma omp parallel for private(numer)
	for(j=0; j < NJ; j++){

		numer = as_scalar(y.row(j)*Exfull - exi.row(j)*Exfull - exi.row(j)*Exa + exixi.row(j)*Exfull - exi.row(j)*Exfull*Epsi(j,0)) + beta_mu/beta_sigma;
		Ebeta(j,0) = numer * Vbeta(j,0);

	}
		
    return;

}



