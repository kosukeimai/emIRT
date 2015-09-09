// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getOnecol_dynIRT(const arma::mat &startlegis,
					const arma::mat &endlegis,
					const int T,
					const int N) {

    arma::mat ones_col(N,T,arma::fill::zeros);
	int i,t;

	for(i=0; i<N; i++){
		for(t=0; t<T; t++){

			if( startlegis(i,0) <= t && endlegis(i,0) >= t){
				ones_col(i,t) = 1.0;
			}
		}
	}

    return(ones_col);

}
