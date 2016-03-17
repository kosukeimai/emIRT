// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
void getEx(arma::mat &Ex,
				   const arma::mat &Vx,
				   const arma::mat &Ealpha,
				   const arma::mat &onesNJ,
                   const arma::mat &exi,
                   const arma::mat &xi,
                   const arma::mat &exixi,
                   const arma::mat &y,
                   const arma::mat &Epsi,
                   const arma::mat &Ebeta,
                   const arma::mat &Vbeta,
                   const arma::mat &i,
                   const double x_mu,
                   const double x_sigma,
                   const int NI,
                   const int NK,
                   const int NJ
                   ) {

	int i_index, j, k;
	double numer;
	arma::mat Ebb(NJ,1);
	arma::mat Ebp(NJ,1);

#pragma omp parallel for
	for(j=0; j<NJ; j++){
		Ebb(j,0) = Ebeta(j,0)*Ebeta(j,0) + Vbeta(j,0);
		Ebp(j,0) = Ebeta(j,0)*Epsi(j,0);
	}

#pragma omp parallel for private(k, numer)
	for(i_index=0; i_index<NI; i_index++){

		numer = 0.0;

		for(k=0; k<NK; k++){
		    if(i(k,0) == i_index) numer += as_scalar( trans(Ebeta)*y.col(k) - trans(Ebeta)*exi.col(k) - trans(Ebp)*exi.col(k)  + trans(Ebeta)*exixi.col(k) - trans(Ebeta)*exi.col(k)*Ealpha(k,0) );
		}

		numer += x_mu/x_sigma;	
		Ex(i_index,0) = numer * Vx(i_index,0);

	}

    return;

}



