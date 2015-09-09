// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEbb_dynIRT(const arma::mat &Eb,
				 const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const int nJ
                 ) {

	int j;
	arma::mat Ebb(nJ,1);

#pragma omp parallel for
	for(j=0; j < nJ; j++){
	    Ebb(j,0) = Eb(j, 0)*Eb(j, 0) + Vb2(1, 1, bill_session(j,0));
	}

    return(Ebb) ;
}
