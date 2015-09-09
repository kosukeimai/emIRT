// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEba_dynIRT(const arma::mat &Ea,
				 const arma::mat &Eb,
				 const arma::cube &Vb2,
                 const arma::mat &bill_session,
                 const int nJ
                 ) {

	int j;
	arma::mat Eba(nJ,1);

#pragma omp parallel for
	for(j=0; j < nJ; j++){
	    Eba(j,0) = Eb(j, 0)*Ea(j, 0) + Vb2(0, 1, bill_session(j,0));
	}

    return(Eba) ;
}
