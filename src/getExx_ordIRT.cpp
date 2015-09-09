// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getExx_ordIRT(const arma::mat &Ex,
                const arma::mat &Vx,
                const int nN
                ){


	int i;
	arma::mat result(nN, 1);

#pragma omp parallel for
    for(i=0; i < nN; i++){

		result(i,0) = Ex(i,0)*Ex(i,0) + Vx(0,0);

    }

    return(result) ;
}
