// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getNlegis_dynIRT(const arma::mat &legis_by_session,
					const int T,
					const int N) {

	int i, counter;
    int maxlegis = legis_by_session.n_cols;
    arma::mat Nlegis(T, 1, arma::fill::zeros) ;
    
    for(i=0; i<T; i++){
		counter = maxlegis;
    	while(legis_by_session(i,counter - 1) == -1){
    		counter = counter - 1;
    	}
		Nlegis(i,0) = counter;

    }
   
    return(Nlegis);
}
