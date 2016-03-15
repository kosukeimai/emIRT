// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getVbeta(arma::mat &Vbeta,
				   const arma::mat &Ex,
				   const arma::mat &Vx,
                   const arma::mat &exi,
                   const arma::mat &i,
				   const double beta_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int j, k;
	double denom;
	arma::mat Exx(NI,1);

#pragma omp parallel for
	for(k=0; k<NI; k++){
		Exx(k,0) = Ex(k,0)*Ex(k,0) + Vx(k,0);
	}

#pragma omp parallel for private(denom, k)
	for(j=0; j<NJ; j++){

		denom = 0.0;
		for(k=0; k<NK; k++){
			denom += exi(j,k)*Exx(i(k,0),0);
		}
		Vbeta(j,0) = 1/(denom + 1/beta_sigma);
	}

    return;

}



