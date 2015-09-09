// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getLast_dynIRT(const arma::mat &bill_session,
					const int T,
					const int J) {

	// Key assumption here and throughout code: bills sorted chronologically before input
	int j, counter;
    arma::mat end_session(T, 1, arma::fill::zeros) ;

	counter=0;    
    for(j=0; j<J; j++){
		if(bill_session(j,0) != counter){
			end_session(counter,0) = j;
			counter = counter + 1;
		}
    }

	end_session(counter,0) = bill_session.n_rows ;

    return(end_session);
}
