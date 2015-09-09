// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getVb_dynIRT(const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const int J
                 ) {

	int j;
	arma::mat Vb(J,1);

#pragma omp parallel for
	for(j=0; j<J; j++){

	    Vb(j,0) = Vb2(1, 1, bill_session(j,0));

	}

    return(Vb) ;
}
