// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getLBS_dynIRT(const arma::mat &startlegis,
					const arma::mat &endlegis,
					const int T,
					const int N) {

    arma::mat legis_by_session(T, N);
	arma::mat LBS;
	int i,j,t,maxlegis;
	maxlegis = 0;
	legis_by_session.fill(-1);

	for(t=0; t<T; t++){
		j=0;
		for(i=0; i<N; i++){

			if( startlegis(i,0) <= t && endlegis(i,0) >= t){
				legis_by_session(t,j) = i;
				j = j + 1;	
			}
		}
		if(j > maxlegis) maxlegis = j;
	}

	LBS = legis_by_session.submat(0,0,T-1,maxlegis-1);
    return(LBS);

}
