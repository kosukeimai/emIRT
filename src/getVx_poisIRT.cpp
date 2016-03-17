// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getVx(arma::mat &Vx,
				   const arma::mat &Ebeta,
				   const arma::mat &Vbeta,
                   const arma::mat &exi,
                   const arma::mat &i,
                   const double x_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int i_index, j, k;
	double denom;
	arma::mat Ebb(NJ,1);

#pragma omp parallel for
	for(j=0; j<NJ; j++){
		Ebb(j,0) = Ebeta(j,0)*Ebeta(j,0) + Vbeta(j,0);
	}

#pragma omp parallel for private(denom, k)
	for(i_index=0; i_index<NI; i_index++){

		denom = 0.0;
		for(k=0; k < NK; k++){
			if(i(k,0) == i_index) denom += as_scalar(trans(Ebb) * exi.col(k));
		}
		Vx(i_index,0) = 1/(denom + 1/x_sigma);

	}

    return;

}



