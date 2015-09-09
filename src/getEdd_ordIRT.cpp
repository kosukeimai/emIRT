// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEdd_ordIRT(const arma::mat &Ex,
                 const arma::mat &Exx,
                 const arma::mat &Eb,
                 const arma::mat &Ebb,
                 const arma::mat &Etau,
                 const arma::mat &Ett,
                 const arma::mat &Ezstar,
                 const arma::mat &Ezzstar,
                 const int N,
                 const int J
                 ) {

	int i,j;
	double alpha, beta;
	arma::mat Edd(J, 1);	

#pragma omp parallel for private(i,j,beta)
	for(j=0; j < J; j++){

		beta = 0.0;

		for(i = 0; i < N; i++){

			// Note: Below does not need modification for either variational or ECM estimate
			// Already made independent when getEtt() and getEbb() updated

	  		beta += Ezzstar(i,j) + Ett(j,0);
//	  		beta += Exx(i,0) * Ebb(j,0) + 2*Ezstar(i,j)*Etau(j,0);
//	  		beta = beta - 2*Ezstar(i,j)*Ex(i,0)*Eb(j,0) - 2*Etau(j,0)*Ex(i,0)*Eb(j,0);
	  		beta += Exx(i,0) * Ebb(j,0) - 2*Ezstar(i,j)*Etau(j,0);
	  		beta = beta - 2*Ezstar(i,j)*Ex(i,0)*Eb(j,0) + 2*Etau(j,0)*Ex(i,0)*Eb(j,0);

		}

		beta = beta/2;
		alpha = 1 + N/2;
		Edd(j,0) = alpha/beta;
	
	}

    return(Edd);

}
