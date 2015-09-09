// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <RcppArmadillo.h>
#include "getEbb_ordIRT.h"

using namespace Rcpp;

// // [[Rcpp::export()]]
arma::mat getEbb_ordIRT(const arma::mat &Eb,
                const arma::mat &Vb,
                const int nJ
                ){


	int j;
	arma::mat result(nJ, 1);

#pragma omp parallel for
    for(j=0; j < nJ; j++){

		// For variational inference estimate
		// result(j,0) = Eb(j,0)*Eb(j,0) + Vb(j,0);

		// For ECM, assume independence here
		result(j,0) = Eb(j,0)*Eb(j,0);


    }

    return(result) ;
}
