// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEx_ordIRT(const arma::mat &Ezstar,
                const arma::mat &Eb,
                const arma::mat &Etau,
                const arma::mat &Vx,
                const arma::mat &Edd,
                const arma::mat &xmu,
                const arma::mat &xsigma,
                const arma::cube &Vb2,
                const int N,
                const int J
                ) {

	int i, j;
//	double tmp;
	arma::mat Ex(N, 1, arma::fill::zeros);
	arma::mat Ebdd(J, 1);
	arma::mat Ebtdd(J, 1);

#pragma omp parallel for
	for(j = 0; j < J; j++){

		//Rcout << Vb2.slice(j)(1,1) << std::endl ;

		Ebdd(j,0) = Eb(j,0)*Edd(j,0); 

		// For variational inference
		// 	Ebtdd(j,0) = ( Eb(j,0) * Etau(j,0) + Vb2.slice(j)(1,0) )*Edd(j,0); 

		// For ECM: E(ba) set as E(b)*E(a)
		Ebtdd(j,0) = ( Eb(j,0) * Etau(j,0)  )*Edd(j,0); 

	}


#pragma omp parallel for private(i,j)
    for(i=0; i < N; i++){

		for(j = 0; j < J; j++){
//	  		Ex(i,0) += Ebdd(j,0)*(Ezstar(i,j) + Etau(j,0));

			// Changed here for joint update.
	  		Ex(i,0) += Ebdd(j,0)*Ezstar(i,j) - Ebtdd(j,0);

		}

		Ex(i,0) = Vx(0,0)*( xmu(0,0)/xsigma(0,0) + Ex(i,0));

    }

    return(Ex) ;

}

